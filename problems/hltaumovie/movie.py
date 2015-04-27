from PIL import Image

file1 = 'movie/orb00160.png'
file2 = '/Users/dtamayo/Desktop/orb00000.png'

hltau = Image.open(file1)
orbits = Image.open(file2)

w,h = hltau.size

ncroptop = 0#125 
ncropleft =0#125 
hltaucrop = hltau.crop((ncropleft,ncroptop,w,h))

w,h = orbits.size
w2,h2 = hltaucrop.size

hltaucrop.thumbnail((w2,h2*h/w), Image.ANTIALIAS)

w,h = 600,600


try:
    hltaucrop.thumbnail((w,h), Image.ANTIALIAS)
    orbits.thumbnail((w,h), Image.ANTIALIAS)
    hltaucrop.save("hltausmall.png", "PNG")
    orbits.save("orbitssmall.png", "PNG")
except IOError:
    print("Resize error")


hltaucrop.putalpha(128)

orbits.paste(hltaucrop, (0,0), hltaucrop)
orbits.save("result.png", "PNG")

import os
os.system("open result.png")

#resize HL Tau image to aspect ratio of output to sidestep non-square images
