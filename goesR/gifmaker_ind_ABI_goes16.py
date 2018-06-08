
import imageio
import numpy as np
nums = [01,02,03,04,05,06,07,010,011,10,11,12,13,14,15,16]
abi_no = []
for n in nums:
    abi_no.append("%02d"%n)


for Band in abi_no:
    with open("/home/sat_ops/goes_r/ind_bands/img_files/img_logfile" + Band + ".txt") as f:
        img_names = f.readlines()
    img_names = [x.strip() for x in img_names] 
    imglen = len(img_names)
    conus_images = []
    midatl_images = []
    dur_vals = []
    
    for i in xrange(1,imglen):
        if i != imglen:
            dur_vals.append(.05)
    dur_vals.append(2)
    
    for i in img_names:
        input_file='/home/sat_ops/goes_r/ind_bands/conus/band' + Band + '/' + str(i)
        conus_images.append(imageio.imread(input_file))
    imageio.mimsave('/home/sat_ops/goes_r/ind_bands/conus/gifs/Band' + Band + '_conus.gif', conus_images, duration=dur_vals)
    
    for i in img_names:
        input_file='/home/sat_ops/goes_r/ind_bands/midAtl/band' + Band + '/' + str(i)
        midatl_images.append(imageio.imread(input_file))
    imageio.mimsave('/home/sat_ops/goes_r/ind_bands/midAtl/gifs/Band' + Band + '_midatlantic.gif', midatl_images, duration=dur_vals)
