import os
import imageio
from PIL import Image

gifs_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots/gifs'
pngs_name = 'inf_1.2_top3'
pngs_dir = os.path.join(gifs_dir, pngs_name)
out_path = os.path.join(gifs_dir, pngs_name + '-512-512-aa-85_top3_jaccard0.gif')

pngs = []
for r, d, f in os.walk(pngs_dir):
    for png in f:
        if '.png' in png:
            pngs.append(os.path.join(r, png))

pngs.sort()

images = []
for png in pngs:
    img = imageio.imread(png)
    img = Image.fromarray(img).resize((512, 512), Image.ANTIALIAS)
    images.append(img)

img.save(out_path, format='GIF', append_images=images, save_all=True, loop=0, quality=85, optimize=True)
