import os
import imageio
from PIL import Image

pngs_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots/dec_paper/maps'
out_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots/dec_paper/maps/gifs'
gif_name = 'inf_1.2_11M_world_large_font_400.gif'
out_path = os.path.join(out_dir, gif_name)

pngs = []
for r, d, f in os.walk(pngs_dir):
    for png in f:
        if '.png' in png:
            pngs.append(os.path.join(r, png))

pngs.sort()

images = []
for png in pngs:
    img = imageio.imread(png)
    img = Image.fromarray(img).resize((400, 400), Image.ANTIALIAS)
    images.append(img)

img.save(out_path, format='GIF', append_images=images, save_all=True, loop=0, quality=85, optimize=True)
