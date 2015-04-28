from PIL import Image


file1 = 'hltau.png'
file2 = "resParticles/orb{0:05d}.png".format(1)

hltau = Image.open(file1)
img = Image.open(file2)

w,h = hltau.size

ncroptop = 160#higher number is down
ncropleft =140#higher is right 
hltaucrop = hltau.crop((ncropleft,ncroptop,w,h))

w,h = img.size
w2,h2 = hltaucrop.size

hltaucrop.thumbnail((w2,h2*h/w), Image.ANTIALIAS)

w,h = 600,600

try:
    hltaucrop.thumbnail((w,h), Image.ANTIALIAS)
    hltaucrop.save("hltausmall.png", "PNG")
except IOError:
    print("Resize error")


import os
os.system("rm ./movie/*")
os.system("rm test.mp4")

FPS = 24
ctr = 0
for i in range(5*FPS):
    file2 = "resParticles/orb{0:05d}.png".format(i)
    img = Image.open(file2)
    img.thumbnail((w,h), Image.ANTIALIAS)
    hltaucrop.putalpha(int(255.-i*255./119.))
    img.paste(hltaucrop, (0,0), hltaucrop)
    img.save("movie/orb{0:05d}.png".format(ctr), "PNG")
    ctr+=1 
for i in range(3*FPS):
    tail = int(i+i**2/8) 
    timg = tail
    base =int(i+i**2/1.33) 
    file1 = 'resParticlesBlurs/orb{0:05d}.png'.format(base)
    img = Image.open(file1)
    img.thumbnail((w,h),Image.ANTIALIAS)
    w,h = img.size
    pixels = img.load()

    for j in range(1,tail):
        file1 = 'resParticlesBlurs/orb{0:05d}.png'.format(base-j)
        img2 = Image.open(file1)
        img2.thumbnail((w,h),Image.ANTIALIAS)
       
        pixels2 = img2.load()

        for j in range(w):
            for k in range(h):
                pixels[j,k] = (max(pixels[j,k][0],pixels2[j,k][0]), max(pixels[j,k][1],pixels2[j,k][1]), max(pixels[j,k][2],pixels2[j,k][2]))

    img.save("movie/orb{0:05d}.png".format(ctr), "PNG")
    ctr += 1
'''
for i in range(8*FPS):
    file1 = 'resOrbits/orb{0:05d}.png'.format(i)
    img = Image.open(file1)
    img.save("movie/orb{0:05d}.png".format(ctr), "PNG")
    ctr += 1
'''    
os.system("ffmpeg -i movie/orb%05d.png -c:v libx264 -pix_fmt yuv420p test.mp4")
os.system("open test.mp4")

