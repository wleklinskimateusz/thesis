# extract nth frame from gif
import os
from PIL import Image


def extract(filename, n):
    im = Image.open(filename)
    try:
        im.seek(n-1)
    except EOFError:
        raise EOFError("There are not %d frames in %s" % (n, filename))
    return im


# if output/iteration/frames doesnt exist, create it
if (not os.path.exists("output/iteration/frames")):
    os.makedirs("output/iteration/frames")

for frame in [1, 4, 6, 9, 13, 24, 54]:
    extract("output/iteration.gif",
            frame).save(f"output/iteration/frames/{frame}.png")
