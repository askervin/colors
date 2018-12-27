#!/usr/bin/env python2
# Copyright 2018 antti.kervinen@gmail.com
# Wavelength-to-RGB conversion is written based on Javascript in
# https://academo.org/demos/wavelength-to-colour-relationship/
"""colors.py - create palette images and recolor images using palette

Usage: colors.py [options]

Options:
  -h, --help                  print help and exit.
  -d, --debug                 increase debug output verbosity.
  -p, --palette PALETTE       use PALETTE.
                              Create a color palette with:
                                c=NUMBER_OF_COLORS
                                l=LIGHT_SHADES_OF_EACH_COLOR
                                d=DARK_SHADES_OF_EACH_COLOR
                                w=LOWEST_WAVELENGTH (default 380)
                                W=HIGHEST_WAVELENGTH (default 760)
                              or a grayscale palette with:
                                g=GRAYSCALES
                              or use colors from existing image:
                                file=image.pnm
                              You can use many palettes by giving many -p.
  -i, --input INPUT.PNM       read image to be recolored from INPUT.PNM
  -o, --output OUTPUT.PNM     write result to OUTPUT.PNM.
                              If also -i is given, OUTPUT.PNM is recolored image.
                              If -i is not given, OUTPUT.PNM consists of
                              colors of all given palettes.

Examples:

Create an image with a 24-color palette:
  $ python colors.py -p 24 -o palette24.pnm

Create an image with 3 palettes:
  1. 32 grayscales
  2. 16 colors with 2 light shades
  3. 16 colors with 2 dark shades.
  $ python colors.py -p g=32 -p c=16,l=2 -p c=16,d=2 -o 3palettes.pnm

Create a blue palette (440 nm...490 nm) with 24 colors and 2 light shades of each color:
  $ python colors.py -p w=440,W=490,c=24,l=2 -o 24x3blues.pnm

Recolor input.pnm using two palettes: 64 shades of gray and colors in 24x3blues.pnm:
  $ python colors.py -i input.pnm -p g=64 -p file=24x3blues.pnm -o result.pnm
"""

import getopt
import math
import re
import sys

g_debug = 0

g_min_wavelen = 380
g_max_wavelen = 760

def output(msg):
    sys.stdout.write(msg)

def error(msg=None, exit_status=1):
    if msg:
        sys.stderr.write('colors.py: %s\n' % (msg,))
    if not exit_status is None:
        sys.exit(exit_status)

def nm_to_rgb(wavelength, gamma=0.8):
    # converts visible light wavelength (nm) to (red, green, blue).
    # red, green and blue are floats between 0..1
    #
    # algorithm converted from javascript to python from:
    # https://academo.org/demos/wavelength-to-colour-relationship/
    factor = 0
    red = 0
    green = 0
    blue = 0
    if (wavelength >= 380) and (wavelength < 440):
        red = - float(wavelength - 440) / (440 - 380)
        green = 0.0
        blue = 1.0
    elif (wavelength >= 440) and (wavelength < 490):
        red = 0.0
        green = float(wavelength - 440) / (490 - 440)
        blue = 1.0
    elif (wavelength >= 490) and (wavelength < 510):
        red = 0.0
        green = 1.0
        blue = - float(wavelength - 510) / (510 - 490)
    elif (wavelength >= 510) and (wavelength < 580):
        red = float(wavelength - 510) / (580 - 510)
        green = 1.0
        blue = 0.0
    elif (wavelength >= 580) and (wavelength < 645):
        red = 1.0
        green = - float(wavelength - 645) / (645 - 580)
        blue = 0.0
    elif (wavelength >= 645) and (wavelength < 781):
        red = 1.0
        green = 0.0
        blue = 0.0
    else:
        red = 0.0
        green = 0.0
        blue = 0.0

    # Let the intensity fall off near the vision limits
    if (wavelength >= 380) and (wavelength < 420):
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif (wavelength >= 420) and (wavelength < 701):
        factor = 1.0
    elif (wavelength >= 701) and (wavelength < 781):
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
    else:
        factor = 0.0

    if (red != 0):
        red = math.pow(red * factor, gamma)

    if (green != 0):
        green = math.pow(green * factor, gamma)

    if (blue != 0):
        blue = math.pow(blue * factor, gamma)

    return (red, green, blue)

def to_pnm(width, height, pixels):
    # convert width x height image that consists of pixels
    # to PNM image format
    max_color = 255
    pnm_header = "P6 %s %s %s\n" % (width, height, max_color)
    pnm_data = []
    for (r, g, b) in pixels:
        pnm_data.append(chr(int(r * max_color)) +
                        chr(int(g * max_color)) +
                        chr(int(b * max_color)))
    return pnm_header + "".join(pnm_data)

def image_to_pnm(image):
    return to_pnm(image['width'], image['height'], image['pixels'])

def palette_of_n_colors(n=12, w_min=380, w_max=760, start_wavelen_offset=0):
    if n == 0:
        return []
    palette = [(0, 0, 0)] * n
    wavelen_change = math.pow(float(w_max)/w_min, 1.0 / n)
    for color_index in xrange(n):
        color_nm = (w_min + start_wavelen_offset) * math.pow(wavelen_change, color_index)
        palette[color_index] = nm_to_rgb(color_nm)
    return palette

def palette_of_n_grayscales(n=12):
    if n == 1:
        return [(0.5, 0.5, 0.5)]
    n_f = float(n-1)
    return [(f/n_f, f/n_f, f/n_f) for f in xrange(n)]

def palette_dark_shades(input_palette, shades=1):
    """returns palette that contains darker shades of colors in input palette"""
    factors = [float(i)/(shades+1.0) for i in xrange(1, shades+1)]
    output_palette = []
    for f in factors:
        for r, g, b in input_palette:
            output_palette.append((f * r, f * g, f * b))
    return output_palette

def palette_light_shades(input_palette, shades=1):
    """returns palette that contains lighter shades of colors in input palette"""
    factors = [float(i)/(shades+1.0) for i in xrange(1, shades+1)]
    output_palette = []
    for f in factors:
        for r, g, b in input_palette:
            output_palette.append((r + f * (1.0 - r),
                                   g + f * (1.0 - g),
                                   b + f * (1.0 - b)))
    return output_palette

def paint_palette(pixels, start_y, height, width, palette):
    if len(palette) == 0:
        return
    palette_width = width / len(palette)
    for y in xrange(start_y, start_y + height):
        x = 0
        for palette_color in palette:
            for palette_x in xrange(palette_width):
                pixels[y*width + x] = palette_color
                x += 1
        while x < width:
            pixels[y*width + x] = palette_color
            x += 1

def create_palette(g=0, c=0, l=0, d=0, w=380, W=760):
    """returns palette based on:
    - g: number of grayscales
    - c: number of colors
    - l: light shades of each color
    - d: dark shades of each color
    - w: minimum wavelength
    - W: maximum wavelength"""

    palette_c = palette_of_n_colors(c, w, W)
    palette_g = palette_of_n_grayscales(g)
    palette_l = palette_light_shades(palette_c, l)
    palette_d = palette_dark_shades(palette_c, d)
    return palette_d + palette_c + palette_l + palette_g

g_palettes = {}

def generate_palette_image(output_pnm, palette_names=None):
    """create image with from given palette or all known palettes"""
    palettes = g_palettes
    width = (g_max_wavelen - g_min_wavelen) * 2
    if palette_names is None:
        palette_names = sorted(g_palettes.keys())
    height = 60 * len(palette_names)
    width = 600
    pixels = [(0,0,0)] * width * height
    y = 0

    for name in palette_names:
        paint_palette(pixels, y, 60, width, palettes[name])
        y += 60

    # full spectrum would be
    #for y in xrange(0, 50):
    #    for x in xrange(width*2):
    #        pixels[y*width + x] = nm_to_rgb(x/2 + g_min_wavelen)
    #        pixels[y*width + x+1] = nm_to_rgb(x/2 + g_min_wavelen)

    open(output_pnm, "wb").write(
        to_pnm(width, height, pixels))

def load_image(image_pnm):
    """returns dictionary {
        'width': width,
        'height': height,
        'pixels': [(Rx0y0, Gx0y0, Bx0y0), (Rx1y0, Gx1y0, Bx1y0), ...]
        } where RGB values are between 0 and 1"""
    file_data = open(image_pnm, "rb").read()
    header_re = re.compile('P6\s([0-9]+)\s([0-9]+)\s([0-9]+)\s')
    m = header_re.match(file_data)
    if not m:
        raise ValueError('file %r is not PNM P6 file' % (image_pnm,))
    width_s, height_s, max_value_s = m.groups()
    if not 0 < int(max_value_s) < 65536:
        raise ValueError('invalid max value in header')
    data = file_data[m.end():]
    width = int(width_s)
    height = int(height_s)
    max_value = float(max_value_s)
    if max_value <= 255.0:
        pixels = [(ord(r) / max_value,
                   ord(g) / max_value,
                   ord(b) / max_value)
                   for r, g, b in zip(*[data[i::3] for i in xrange(3)])]
    else:
        rr = zip(data[0::6], data[1::6])
        gg = zip(data[2::6], data[3::6])
        bb = zip(data[4::6], data[5::6])
        pixels = [((ord(rr[i][0])*256 + ord(rr[i][1])) / max_value,
                   (ord(gg[i][0])*256 + ord(gg[i][1])) / max_value,
                   (ord(bb[i][0])*256 + ord(bb[i][1])) / max_value)
                   for i in xrange(len(bb))]

    if len(pixels) != width * height:
        raise ValueError('invalid pixel array size, got %s pixels, expected %s x %s = %s' %
                         (len(pixels), width, height, width*height))
    return {'width': width,
            'height': height,
            'pixels': pixels}

def load_palette_file(image_pnm):
    """returns palette containing colors in input image"""
    image = load_image(image_pnm)
    found_colors = set()
    for pixel in image['pixels']:
        if not pixel in found_colors:
            found_colors.add(pixel)
    return sorted(found_colors)

def recolor_image(input_pnm, palette, output_pnm):
    """reads input_pnm, recolors it using palette, saves output to output_pnm"""
    image = load_image(input_pnm)
    color_hash = {}
    pixels = image['pixels']
    for pixel_index in xrange(len(pixels)):
        pixel_color = pixels[pixel_index]
        try:
            best_pal_color = color_hash[pixel_color]
        except KeyError:
            r, g, b = pixel_color
            min_error = 4.0 # definitely bigger than 1.0**2 + 1.0**2 + 1.0**2
            best_pal_color = (0.0, 0.0, 0.0)
            palette_error = []
            for (R, G, B) in palette:
                pal_error = (r-R)**2 + (g-G)**2 + (b-B)**2
                if pal_error < min_error:
                    min_error = pal_error
                    best_pal_color = (R, G, B)
            color_hash[pixel_color] = best_pal_color
        pixels[pixel_index] = best_pal_color
    open(output_pnm, "wb").write(
        image_to_pnm(image))

def main(argv):
    global g_debug
    opt_input = None
    opt_output = None
    opt_palette = []
    opts, remainder = getopt.getopt(
        argv,
        "hdi:p:o:",
        ["debug", "help",
         "input=", "palette=", "output="])
    for opt, arg in opts:
        if opt in ["-h", "--help"]:
            output(__doc__)
            error(exit_status=0)
        elif opt in ["-d", "--debug"]:
            g_debug += 1
        elif opt in ["-i", "--input"]:
            opt_input = arg
        elif opt in ["-o", "--output"]:
            opt_output = arg
        elif opt in ["-p", "--palette"]:
            if "file=" in arg:
                g_palettes[arg] = load_palette_file(arg.split("file=")[-1])
            elif "=" in arg:
                specs = {}
                try:
                    for key_eq_val in arg.split(","):
                        key, val = key_eq_val.split("=")
                        if not key in "gcldwW":
                            raise ValueError()
                        specs[key] = int(val)
                except Exception:
                    error('invalid --palette %r, use syntax g=GRAYSCALES,c=COLORS,l=LIGHT_COLOR_SHADES,d=DARK_COLOR_SHADES,w=MIN_NM,W=MAX_NM' % (arg,))
                g_palettes[arg] = create_palette(**specs)
            elif not arg in g_palettes:
                error('invalid palette name %r, see --help' % (arg,))
            opt_palette.append(arg)

    if opt_input is None and opt_palette == [] and not opt_output is None:
        generate_palette_image(opt_output)
        error(exit_status=0)
    elif opt_input is None and (opt_palette) and (not opt_output is None):
        generate_palette_image(opt_output, palette_names=opt_palette)
        error(exit_status=0)
    elif opt_input is None:
        error('missing --input INPUT.PNM')
    elif opt_palette is None:
        error('missing --palette NAME')
    elif opt_output is None:
        error('missing --output OUTPUT.PNM')

    combined_palette = []
    for palette_name in opt_palette:
        combined_palette.extend(g_palettes[palette_name])
    recolor_image(opt_input, list(set(combined_palette)), opt_output)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Exception, e:
        error(str(e) + " (%s)" % (type(e),))
