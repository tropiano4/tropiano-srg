# Scratch Work

# Test codes in progress in this script.
# Last thing tested: Re-write Wendt potential/omega files to new kvnn numbers


from os import chdir, getcwd, rename, walk

cwd = getcwd()

#kvnn_new = '900'
#kvnn_new = '901'
kvnn_new = '902'

chdir('Potentials/vsrg_macos/vsrg_kvnn_%s_lam12.0_kmax30_kmid4_ntot120' % \
      kvnn_new)

for root, dirs, files in walk('.'):
    
    for name in files:
        
        name_list = name.split('_')
        
        if name_list[0] in ['omega', 'vnn', 'vsrg']:
            
            name_list[3] = kvnn_new
            new_name = '_'.join(name_list)
            
            rename(name, new_name) 
            
chdir(cwd)