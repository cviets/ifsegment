from .modules.cyto_segment import cyto_segment_folder

def main(input, output, channel):
    return cyto_segment_folder(input, output, channel)