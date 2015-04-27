from PIL import Image

w,h = 600,600

file1 = 'resParticles/orb{0:05d}.png'.format(0)
img = Image.open(file1)
img.thumbnail((w,h),Image.ANTIALIAS)
w,h = img.size
pixels = img.load()

for i in range(1,3):
    file1 = 'resParticles/orb{0:05d}.png'.format(i)
    img2 = Image.open(file1)
    img2.thumbnail((w,h),Image.ANTIALIAS)
   
    pixels2 = img2.load()
    for j in range(w):
        for k in range(h):
            old = pixels[j,k]
            pixels[j,k] = (max(pixels[j,k][0],pixels2[j,k][0]), max(pixels[j,k][1],pixels2[j,k][1]), max(pixels[j,k][2],pixels2[j,k][2]))
            if pixels[j,k] != old:
                print("changed")
img.save("result.png", "PNG")

import os
os.system("open result.png")

#resize HL Tau image to aspect ratio of output to sidestep non-square images
