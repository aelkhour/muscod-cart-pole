# \file   default.plot
# \author Christian Kirches
# \author $LastChangedBy: ckirches $
# \date   July 20, 2007
# \date   $LastChangedDate: 2007-12-14 14:49:35 +0100 (Fri, 14 Dec 2007) $
#
# $Id: default.plot.template 565 2007-12-14 13:49:35Z ckirches $
#
# \brief  Default PLOT configuration for MUSCOD-II
#
#         Everything herein is case-sensitive !
#

# General configuration

leave_windows_open       off       # Leave windows open when multiple runs are done (MS-MINTOC) ?
draw_shooting_nodes      on        # draw multiple-shooting nodes in xd and xa trajectories
draw_continuous_controls off       # connect discontinuous controls across shooting intervals
draw_bounds              on        # draw lower/upper bounds in all panels

# Screen output configuration

screen_enabled	   on        # enable output to screen ?
screen_width       15        # inch = pixels / dpi (10.6 for 1024 pixels, 13.3 for 1280 pixels, 16.6 for 1600 pixels)
screen_aspect      0.625      # height/width (0.75 for a 4:3 screen, 0.625 for a 16:10 widescreen)
screen_hostname    ""        # PGPLOT x server name, e.g., jim.iwr.uni-heidelberg.de
screen_display     0         # Usually 0
screen_screen      0         # Usually 0 

# Postscript output configuration

ps_enabled	   	   off      # enable output to postscript file ? 
ps_append_file	   on        # All iterations into the same file (multiple pages) ?
ps_width           7.5       # A4 paper width in inches (7.5 for portrait, 10.6 for landscape)
ps_aspect          1.4142    # height/width (1.4142 for portrait, 0.7071 for landscape)
ps_color           on        # color or grayscale postscript ?

# Gif output configuration

gif_enabled        off       # enable output to gif file ?
gif_width          8.0       # width in inches
gif_aspect         0.75      # height in inches

# gif_append_file  is not supported by the PGPLOT GIF driver
# gif_color        is not supported by the PGPLOT GIF driver

# end of file
