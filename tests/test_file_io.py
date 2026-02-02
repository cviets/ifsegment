from src.ifsegment.modules.io_utils import read_czi, get_well_from_file, get_czi_in_folder
import numpy as np

def test_read_czi(inp, save_to):
    img = read_czi(inp)
    print(type(img), img.shape)
    np.save(save_to, img)

def test_well_from_file(filename):
    print(get_well_from_file(filename))

def test_get_czi_from_folder(path_to_czi_files):
    return get_czi_in_folder(path_to_czi_files)

def main():
    path = "/media/cviets/Chris2/Exp008J-01"
    return test_get_czi_from_folder(path)

if __name__ == '__main__':
    print(main())