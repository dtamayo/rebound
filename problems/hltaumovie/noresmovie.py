from PIL import Image


file1 = 'hltau.png'
file2 = "noresOrbitsFade/orb{0:05d}.png".format(0)

hltau = Image.open(file1)
img = Image.open(file2)

w,h = hltau.size

ncroptop = 0#higher number is down
ncropleft =0#higher is right 
hltaucrop = hltau.crop((ncropleft,ncroptop,w,h))

w,h = img.size
w2,h2 = hltaucrop.size

hltaucrop.thumbnail((w2,h2*h/w), Image.ANTIALIAS)

w,h = 550,550

try:
    hltaucrop.thumbnail((w,h), Image.ANTIALIAS)
    hltaucrop.save("hltausmall.png", "PNG")
except IOError:
    print("Resize error")


import os
os.system("rm ./noresmovie/*")
os.system("rm nores.mp4")

FPS = 24
ctr = 0

for i in range(5*FPS):
    file2 = "noresOrbitsFade/orb{0:05d}.png".format(i)
    img = Image.open(file2)
    img.thumbnail((w,h), Image.ANTIALIAS)
    hltaucrop.putalpha(int(160.-i*160./119.))
    img.paste(hltaucrop, (0,0), hltaucrop)
    img.save("noresmovie/orb{0:05d}.png".format(ctr), "PNG")
    ctr+=1

ctr1 = ctr
for i in range(5*FPS):
    file1 = 'noresOrbitsFast/orb{0:05d}.png'.format(ctr-ctr1)
    img = Image.open(file1)
    img.thumbnail((w,h), Image.ANTIALIAS)
    img.save("noresmovie/orb{0:05d}.png".format(ctr), "PNG")
    ctr+=1

ctr2 = ctr
for i in range(1*FPS):
    file1 = 'noresOrbitsZoom/orb{0:05d}.png'.format(ctr-ctr2)

    img = Image.open(file1)
    img.thumbnail((w,h),Image.ANTIALIAS)
    w,h = img.size

    file1 = 'noresEject/orb{0:05d}.png'.format(ctr-ctr2)
    img2 = Image.open(file1)
    img2.thumbnail((w,h),Image.ANTIALIAS)
       
    img.putalpha(int(255.-i*255./23.))
    img2.paste(img, (0,0), img)
    img2.save("noresmovie/orb{0:05d}.png".format(ctr), "PNG")
    ctr += 1

for i in range(7*FPS):
    file1 = 'noresEject/orb{0:05d}.png'.format(i+24)
    img = Image.open(file1)
    img.save("noresmovie/orb{0:05d}.png".format(ctr), "PNG")
    ctr += 1

os.system("ffmpeg -i noresmovie/orb%05d.png -c:v libx264 -pix_fmt yuv420p nores.mp4")
os.system("open nores.mp4")
