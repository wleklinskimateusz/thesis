# write a script that changes all the filenames in a given directory

import os

DIR = "output/anim/frames"


def change_name():
    files = os.listdir(DIR)
    for file in files:
        os.rename(
            os.path.join(DIR, file),
            os.path.join(
                DIR, file
                .replace(".png", "PLACEHOLDER")
                .replace(".", "_")
                .replace("PLACEHOLDER", ".png")
            )
        )


change_name()
