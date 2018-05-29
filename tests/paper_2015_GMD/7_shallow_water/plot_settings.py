import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def ticks_changes(ax):
    # removing some ticks' labels
    for i, tick in enumerate(ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks()):
        if i % 2 != 0:
            tick.label1On = False

    # changing ticks' size
    for item in plt.xticks()[1] + plt.yticks()[1]:
        item.set_fontsize(15)
