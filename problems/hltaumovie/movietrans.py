from PIL import Image
import os
os.system("rm -rf movietest/*")
os.system("rm -rf movie/*")
os.system("rm test.mp4")

FPS = 24
for i in range(7*FPS,8*FPS):
    file1 = 'resPtrans/orb{0:05d}.png'.format(i)
    file2 = 'resOrbitsTrans/orb{0:05d}.png'.format(i-7*FPS)
    img = Image.open(file1)
    img2 = Image.open(file2)
    img2.paste(img, (0,0), img)
    img2.save("movietest/orb{0:05d}.png".format(i),"PNG")

os.system("cp moviebackup/* movie/")
os.system("cp movietest/* movie/")
os.system("ffmpeg -i movie/orb%05d.png -c:v libx264 -pix_fmt yuv420p test.mp4")
os.system("open test.mp4")
