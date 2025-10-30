from .modules.nucleus_segment import segment_folder

def main(input, output, channel_n, channel_c, mode):
    return segment_folder(input, output, channel_n, channel_c, mode)