import numpy as np
from typing import Tuple, List
from .io_utils import get_czi_in_folder, get_well_from_file, get_mask_path, read_czi, read_tiff, write_to_csv
from .normalizations import zstack
from tqdm import tqdm

def quantify_channels(image: np.ndarray[np.float64], mask: np.ndarray[np.float64]) -> Tuple[np.ndarray[np.bool_]]:
    """
    Quantify nuclear and cytoplasmic protein fluorescence

    Parameters 
    -----------
    image : np.ndarray[np.float64]
        Input image as numpy array (CYX with C containing channels of interest)
    mask : np.ndarray[np.float64]
        Trinary mask (2D array, YX, 0 = background, 1 = nucleus, 2 = cytoplasm)

    Returns 
    -----------
    fluor : np.ndarray[np.float64]
        C-by-2 array containing [nuclear, cytoplasmic] fluorescence for each protein
    """
    # NOTE: it is important not to modify the pixel values in the image itself since we don't want to mess with the data
    nuc_mask = mask==1
    cyto_mask = mask==2

    if image.ndim == 2:
        image = np.expand_dims(image, 0)

    fluor = np.zeros(shape=(len(image), 3))

    # compute mean intensities for each channel
    for i, channel in enumerate(image):
        fluor[i, 0] = np.mean(channel[np.logical_or(nuc_mask, cyto_mask)])
        fluor[i, 1] = np.mean(channel[nuc_mask])
        fluor[i, 2] = np.mean(channel[cyto_mask])

    return fluor

def quantify_folder(path_to_images: str, path_to_masks: str, path_to_save: str, channels: List[int], mode:str):
    """
    Parameters 
    -----------
    path_to_images : str
        Path to czi images
    path_to_masks : str
        Path to tiff trinary masks
    path_to_save : str
        Path to save data (CSV file)
    channels : List[int]
        List of channel indices (0-indexed) to quantify fluorescence
    mode : str
        Mode for z-projecting (max or mean)
    
    Returns 
    -----------
    None
    """
    image_paths = get_czi_in_folder(path_to_images)

    # for each channel, we measure total, nuclear, cytoplasmic, and N/C ratio
    data = np.zeros(shape=(len(image_paths), 1+4*len(channels)),dtype=object)
    
    for idx, image_path in enumerate(tqdm(image_paths, desc="Measuring all images")):
        well_name = get_well_from_file(image_path)
        mask_path = get_mask_path(path_to_masks, well_name)
        image_czi = read_czi(image_path)
        image_dask = image_czi.get_image_dask_data("CZYX", C=channels)
        image_numpy = image_dask.compute()
        image = zstack(image_numpy, 1, mode)
        mask = read_tiff(mask_path)
        data[idx, 0] = well_name

        fluor = quantify_channels(image, mask)
        for idx_j in range(len(channels)):
            total = fluor[idx_j, 0]
            nuclear = fluor[idx_j, 1]
            cytoplasmic = fluor[idx_j, 2]
            data[idx, 4*idx_j+1] = total
            data[idx, 4*idx_j+2] = nuclear
            data[idx, 4*idx_j+3] = cytoplasmic
            data[idx, 4*idx_j+4] = nuclear/cytoplasmic

    header = ["Well"]
    for channel in channels:
        channel_string = "Ch"+str(channel)
        header += [channel_string+"_TOT", channel_string+"_N", channel_string+"_C", channel_string+"_N/C"]
    data = np.vstack((header, data))
    write_to_csv(data, path_to_save)