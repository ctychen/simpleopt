import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from stl import mesh
import imageio
import ffmpeg


def stl_to_image(stl_file, image_file):
    print(f"Processing {stl_file}...")
    stl_mesh = mesh.Mesh.from_file(stl_file)

    figure = plt.figure()
    axes = Axes3D(figure)

    x = stl_mesh.vectors[:,:,0]
    y = stl_mesh.vectors[:,:,1]
    z = stl_mesh.vectors[:,:,2]

    axes.plot_trisurf(x.flatten(), y.flatten(), z.flatten())
    axes.set_box_aspect([np.ptp(a) for a in (x, y, z)])

    print(f"Saving image {image_file}...")
    plt.savefig(image_file)
    print(f"Image saved {image_file}.")


def generate_video(image_folder, output_video):
    print("Generating video...")
    (
        ffmpeg
        .input(os.path.join(image_folder, '*.png'), pattern_type='glob', framerate=25)
        .output(output_video)
        .run()
    )
    print(f"Video saved as {output_video}.")


def main(stl_folder, output_video):
    image_folder = 'images'
    os.makedirs(image_folder, exist_ok=True)

    stl_files = glob.glob(os.path.join(stl_folder, '*.stl'))

    for i, stl_file in enumerate(stl_files):
        image_file = os.path.join(image_folder, f'image{i}.png')
        try:
            stl_to_image(stl_file, image_file)
        except Exception as e:
            print(f"Failed to process {stl_file}. Reason: {e}")
        else:
            plt.close()  # Close the figure only if the image was successfully saved

    generate_video(image_folder, output_video)


if __name__ == '__main__':
    stlFolder = '/Users/cchen/Desktop/simpleopt/src/pyramidtest013'
    main(stlFolder, 'output.mp4')
