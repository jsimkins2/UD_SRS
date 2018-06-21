
with open("/home/sat_ops/goes_r/lightning/img_logfile.txt") as f:
    img_names = f.readlines()
img_names = [x.strip() for x in img_names] 

imglen = len(img_names)


import imageio
import numpy as np
images = []
images2 = []

dur_vals = []
for i in xrange(0,imglen - 1):
    dur_vals.append(.07)
    
dur_vals.append(1.5)
#print dur_vals
new_seq = range(0,imglen)
from collections import OrderedDict
#new_seq = sorted(new_seq, key=int, reverse=True)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/lightning/multi_con/' + img_names[i]
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/lightning/lightning_conus.gif', images, duration=dur_vals)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/lightning/multi_mid/' + img_names[i]
    images2.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/lightning/lightning_midatlantic.gif', images2, duration=dur_vals)




