# Looks at the sample def and determines which feature-counter calls can be disabled
import os
import martian
__MRO__ = """
stage DISABLE_SPATIAL_STAGES(
    in  path     spot_image_path,
    in  path     tissue_image_path,
    out bool     disable_spatial,
    src py     "stages/spatial/disable_spatial",
)
"""

def main(args, outs):
    if not args.tissue_image_path:      # spot image path will become optional
        outs.disable_spatial = True
    else:
        if args.spot_image_path is not None:
            if not os.path.exists(args.spot_image_path):
                martian.exit("Spot image file not found: {}".format(args.spot_image_path))
        if not os.path.exists(args.tissue_image_path):
            martian.exit("Tissue image not found: {}".format(args.tissue_image_path))
        outs.disable_spatial = False
