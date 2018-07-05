# conus nexrad gifmaker

with open("/home/sat_ops/goesR/radar/imgconus_logfile.txt") as f:
    img_names = f.readlines()
img_names = [x.strip() for x in img_names] 

imglen = len(img_names)

import imageio
import numpy as np
images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.07)
dur_vals.append(2)

for i in img_names:
    input_file='/home/sat_ops/goesR/radar/imgconus/' + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goesR/radar/nexrad_conus.gif', images, duration=dur_vals)


