from .cyto_segment import cyto_segment_czi
from aicsimageio import AICSImage
import numpy as np
from .normalizations import minmax_percentile, clip_512, zstack
from .io_utils import get_czi_in_folder, read_czi, get_well_from_file, write_tiff
from tqdm import tqdm
from skimage.measure import label, regionprops
from skimage.morphology import remove_small_objects, remove_small_holes, diamond
from skimage.segmentation import find_boundaries
from typing import Tuple
import copy
from scipy.ndimage import distance_transform_edt
from scipy.ndimage import binary_closing

def nucleus_preprocess(czi_image: AICSImage, nuclear_channel: int, cyto_channel: int, mode:str) -> np.ndarray[np.bool_]:
    """
    Generates nuclear mask from input AICS Image (from read_czi)

    Parameters 
    -----------
    czi_image : AICSImage
        AICSImage object containing CZI image information
    nuclear_channel : int
        Nuclear marker channel number (0-indexed) (eg, Hoechst)
    cyto_channel : int
        Cytoplasmic marker channel number (0-indexed), used to mask nuclei with associated cytoplasm

    Returns
    ---------
    masks : np.ndarray[np.bool_]
        Segmentation array containing cytoplasmic and nuclear masks
    """

    cyto_mask = cyto_segment_czi(czi_image, cyto_channel, mode)

    img_dask = czi_image.get_image_dask_data("ZYX", T=0, C=nuclear_channel)
    img_numpy = img_dask.compute()
    img_numpy = zstack(img_numpy, axis=0, mode=mode)

    # HARD CODED: normalize 0 percentile to -1 and 98 percentile to 1
    img_numpy = clip_512(img_numpy)
    img_numpy = minmax_percentile(img_numpy, 0, 98)

    return img_numpy, cyto_mask

def nuc_segment_array(
        image: np.ndarray[np.float64], 
        cyto_mask: np.ndarray[np.bool_]
        ) -> Tuple[np.ndarray[np.bool_], int]:
    """
    Create nuclear segmentation of input 2D image 

    Parameters 
    -----------
    image: np.ndarray[np.float64]
        Input image as numpy array
    cyto_mask: np.ndarray[np.bool_]
        Cytoplasm mask (used to determine valid nuclei)
    
    Returns 
    -----------
    mask_out : np.ndarray[np.bool_]
        Mask of image
    num_cells: int
        Number of nuclei counted
    """

    assert image.ndim == 2

    # we are interested in the non-one values 
    # (want to ignore super outlier-y bright spots)
    mu = np.mean(image[image < 1])
    std = np.std(image[image < 1])

    # hard coded: mask according to mean +/- 5 std for non-bright spots
    thresh = min(0, mu + 5*std)
    mask = image > thresh

    mask = remove_small_objects(mask, min_size=300, connectivity=2)
    # mask = remove_small_holes(mask, 1200, connectivity=2)

    # remove nuclei not touching cell cytoplasm (usually means no neurites)
    cyto_dt = distance_transform_edt(~cyto_mask)
    mask_out = np.zeros_like(mask, dtype=np.bool_)

    labeled_nuclei = label(mask)
    regions = regionprops(labeled_nuclei)
    num_cells = 0
    for region in tqdm(regions, desc="Validating nuclei"):

        cur_nucleus = region.image
        boundary = find_boundaries(cur_nucleus, mode="outer")
        temp_dt = cyto_dt[region.slice]
        _temp_dt = copy.deepcopy(temp_dt)
        temp_dt[~boundary] = -1
        distances = temp_dt[temp_dt != -1]

        # majority rules
        if len(distances[distances==0]) >= 0.5 * len(distances):
            num_cells += 1
            mask_out[region.slice] = cur_nucleus

        cyto_dt[region.slice] = _temp_dt
    
    mask_out = binary_closing(mask_out, diamond(1), iterations=1)
    mask_out = remove_small_holes(mask_out, 100, connectivity=2)
    return mask_out, num_cells

def remove_unconnected_cyto(cyto_mask:np.ndarray[np.bool_], nuc_mask: np.ndarray[np.bool_]) -> np.ndarray[np.bool_]:
    """
    Remove cytoplasmic objects not touching or overlapping nuclear objects
    """
    mask_out = np.zeros_like(cyto_mask)
    labeled_cyto = label(cyto_mask)
    regions = regionprops(labeled_cyto)
    nuc_dt = distance_transform_edt(~nuc_mask)
    for region in tqdm(regions, desc="Validating cytoplasm"):

        cur_cyto = region.image
        temp_dt = nuc_dt[region.slice]
        _temp_dt = copy.deepcopy(temp_dt)
        temp_dt[~cur_cyto] = -1
        distances = temp_dt[temp_dt != -1]

        if min(distances) <= 1:
            mask_out[region.slice] = cur_cyto

        nuc_dt[region.slice] = _temp_dt

    return mask_out

def fill_holes_trinary(trinary_mask):
    binary_mask = trinary_mask > 0
    binary_mask_filled = remove_small_holes(binary_mask, 200)
    holes = np.logical_and(binary_mask_filled, ~binary_mask)
    mask_out = trinary_mask
    mask_out[holes] = 2
    return mask_out

def segment_folder(path_to_folder: str, output_folder: str, nuclear_channel: int, cyto_channel: int, mode:str) -> None:
    czi_paths = get_czi_in_folder(path_to_folder)
    cell_counts = [None]*len(czi_paths)
    for i, czi_path in enumerate(tqdm(czi_paths, desc="Masking cells")):
        img = read_czi(czi_path)
        well_name = get_well_from_file(czi_path)
        preprocessed_img, cyto_mask = nucleus_preprocess(img, nuclear_channel, cyto_channel, mode)
        write_tiff(preprocessed_img, output_folder, well_name + "_PREPROCESSED_NUC")
        nuc_mask, num_cells = nuc_segment_array(preprocessed_img, cyto_mask)
        cyto_mask = remove_unconnected_cyto(cyto_mask, nuc_mask)
        mask = np.zeros_like(nuc_mask, dtype=np.float64)

        # trinary image: nuc = 1, cyto = 2
        mask[cyto_mask] = 2
        # nuc_mask should override cyto_mask
        mask[nuc_mask] = 1

        mask = fill_holes_trinary(mask)
        
        write_tiff(mask, output_folder, well_name)
        cell_counts[i] = num_cells

    return cell_counts
        