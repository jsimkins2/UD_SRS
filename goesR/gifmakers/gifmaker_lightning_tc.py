
with open("/home/sat_ops/goes_r/lightning/add_clouds/img_logfile.txt") as f:
    img_names = f.readlines()
img_names = [x.strip() for x in img_names] 

imglen = len(img_names)


import imageio
import numpy as np
images3 = []
images4 = []

dur_vals = []
for i in xrange(0,imglen - 1):
    dur_vals.append(.07)
    
dur_vals.append(1.5)
#print dur_vals
new_seq = range(0,imglen)
from collections import OrderedDict
#new_seq = sorted(new_seq, key=int, reverse=True)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/lightning/add_clouds/im_tc_con/' + img_names[i]
    images3.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/lightning/add_clouds/lightning_true_color_conus.gif', images3, duration=dur_vals)

for i in new_seq:
    input_file='/home/sat_ops/goes_r/lightning/add_clouds/im_tc_mid/' + img_names[i]
    images4.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/lightning/add_clouds/lightning_true_color_midatlantic.gif', images4, duration=dur_vals)
