import numpy as np
from astropy.io import fits as pyfits
import argparse


def powerlaw(x, index, norm):
    print(x,index, norm)
    return norm * x**(-1. * index)

parser=argparse.ArgumentParser()

parser.add_argument('--rmf', dest="rmf", type=str, default="MSFCRMF_Stage2.rmf",help="The RMF used to estimate the count rate")
parser.add_argument('--arf', dest='arf', type=str, default="artxc_arf_2048.arf",help="The ARF used to estimate the count rate")
parser.add_argument('--elo', dest='elo', type=float, default=4.0,help="Minumum energy to include in count rate")
parser.add_argument('--ehi',dest='ehi',type=float,default=15.0,help="Maximum energy to include in count rate")




args = parser.parse_args()



rmffile = args.rmf
arffile = args.arf
elo = args.elo
ehi = args.ehi

rmf_hdulist = pyfits.open(rmffile)

rmf_primary = rmf_hdulist[0]



channel_hdu = rmf_hdulist[2]  #This hdu defines the channels and bounds to divide up the intrinsic input spectrum


output_channel_arr = channel_hdu.data['channel']
output_channel_emin = channel_hdu.data['e_min']
output_channel_emax = channel_hdu.data['e_max']

output_channel_emid = (output_channel_emax + output_channel_emin)/2.


rmf_hdu = rmf_hdulist[1]

rmf_data = rmf_hdu.data

rmf_header = rmf_hdu.header

n_channels = rmf_header['naxis2']  #This is the number of input energy channels


#print(n_channels)
rmf_matrix = rmf_data['MATRIX']  #This matrix is n_channels x max(allchannel_arr), and transforms an intrinsic spectrum into a redistributed detector spectrum 

matrix_list = []
#Put it into the correct numpy matrix format for dot products
for j in range(n_channels):

    this_energ_lo = rmf_data['ENERG_LO'][j]
    this_energ_hi = rmf_data['ENERG_HI'][j]
    this_row = rmf_matrix[j]

    matrix_list.append(this_row)
    
matrix_arr = np.array(matrix_list)

arf_hdulist = pyfits.open(arffile)
arf_hdu =arf_hdulist[1]
#print(arf_hdulist, arf_hdu, arf_hdu.data)
arf_data = arf_hdu.data
arf_elo = arf_data['energ_lo']
arf_eup = arf_data['energ_hi']

arf_emid = 0.5*(arf_elo + arf_eup)
arf_de = 0.5*(arf_eup - arf_elo)
arf_area = arf_data['specresp']


#print(np.shape(matrix_arr))

input_channel_emin = rmf_data['ENERG_LO']
input_channel_emax = rmf_data['ENERG_HI']

input_channel_emid = (input_channel_emax + input_channel_emin)/2.

model_dict = {"powerlaw" : ['gamma','norm'] }
parameter_value_list = []
parameter_value_dict = {}
print("Choose a model to simulate the count rate for. Your choices (and the parameters you need to input after that) are listed below")
for this_key in model_dict.keys():
    print("Model: ", this_key, "|| Model Parameters: ", model_dict[this_key])
model_to_use = input()

if not(model_to_use in model_dict.keys()):
    print("%s is not a valid model. Goodbye")
    quit()

for this_parameter in model_dict[model_to_use]:
    print("We now need to provide inputs for the model parameters")
    this_par_value = input("Parameter %s " %(this_parameter))
    this_par_value = float(this_par_value)
    parameter_value_list.append(this_par_value)
    parameter_value_dict[this_parameter] = this_par_value
    
model_flux = eval(model_to_use)(arf_emid, *parameter_value_list) 

energy_mask_lo = arf_emid >= args.elo
energy_mask_up = arf_emid <= args.ehi

energy_mask = np.logical_and(energy_mask_lo,energy_mask_up)

#arf_de has the half width of the energy bin, hence the factor of 2. 
#1.60218e-9 is the conversion factor between keV and erg
flux_model = np.sum(model_flux[energy_mask]* arf_emid[energy_mask] * 2.*arf_de[energy_mask]) * 1.60218e-9
print("Flux of model: %.3e in the %.1f - %.1f keV energy band" %(flux_model, args.elo, args.ehi))
print(args.elo, args.ehi)
#for(ee,ff, aa,dd) in zip(arf_emid[energy_mask],model_flux[energy_mask], arf_area[energy_mask], arf_de[energy_mask]) : print(ee,ff,aa,dd)

cts_input_flux = model_flux * arf_area * 2.*arf_de

cts_output_flux = np.dot(matrix_arr.transpose(), cts_input_flux)

energy_mask_out_lo = output_channel_emid >= args.elo
energy_mask_out_up = output_channel_emid <= args.ehi

energy_mask_out = np.logical_and(energy_mask_out_lo,energy_mask_out_up)

ctrate = np.sum(cts_output_flux[energy_mask_out])
print("Model Count Rate: %.3e in the %.1f-%.1f keV energy band" %(ctrate, args.elo, args.ehi))
