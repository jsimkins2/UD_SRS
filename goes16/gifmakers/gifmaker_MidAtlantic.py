with open("/home/sat_ops/goes_r/cloud_prod/noaa_format/MAimg_logfile.txt") as f:
    img_names = f.readlines()
img_names = [x.strip() for x in img_names] 

imglen = len(img_names)

import imageio
import numpy as np
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.1)
dur_vals.append(2)
#print dur_vals

for i in img_names:
    input_file='/home/sat_ops/goes_r/cloud_prod/noaa_format/image_midatlantic/' + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/cloud_prod/noaa_format/true_color_midatlantic.gif', images, duration=dur_vals)