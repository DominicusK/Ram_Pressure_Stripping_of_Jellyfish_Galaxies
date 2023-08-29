#Find the line width of your Latex document in pixels. This can be done, temporarily inserting somewhere in your document, \showthe\columnwidth
def set_size(width, fraction=1, subplots=(2, 2)):

    width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)
#
figsize=set_size(width=508, fraction=1, subplots=(2, 2))
print(figsize)
