import sys
import imageio
import numpy as np

num = int(sys.argv[1])

if num==1:
    site = 'KDOX'
    image_dir = 'image_kdox_goes/'
    data_dir = 'data_kdox/'

if num==2:
    site = 'KDIX'
    image_dir = 'image_kdix_goes/'
    data_dir = 'data_kdix/'

if num==3:
    site = 'KLWX'
    image_dir = 'image_klwx_goes/'
    data_dir = 'data_klwx/'

with open("/home/sat_ops/goes_r/nexrad/imggoesnxrd_logfile.txt") as f:
    img_names = f.readlines()
img_names = [x.strip() for x in img_names] 

imglen = len(img_names)


images = []
dur_vals = []
for i in xrange(1,imglen):
    if i != imglen:
        dur_vals.append(.1)
dur_vals.append(2)
#print dur_vals

for i in img_names:
    input_file='/home/sat_ops/goes_r/nexrad/' + image_dir + str(i)
    images.append(imageio.imread(input_file))
imageio.mimsave('/home/sat_ops/goes_r/nexrad/kdox_nexrad_goes16.gif', images, duration=dur_vals)
