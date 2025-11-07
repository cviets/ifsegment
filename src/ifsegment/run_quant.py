from .modules.protein_quantification import quantify_folder

def main(input, masks, output, channels, mode):
    if isinstance(channels, int):
        channels = [channels]
    quantify_folder(input, masks, output, channels, mode)