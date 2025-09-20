from src.ifsegment.modules.io_utils import read_czi, get_well_from_file
import numpy as np

def test_read_czi(inp, save_to):
    img = read_czi(inp)
    print(type(img), img.shape)
    np.save(save_to, img)

def test_well_from_file(filename):
    print(get_well_from_file(filename))

def main():
    return

if __name__ == '__main__':
    main()