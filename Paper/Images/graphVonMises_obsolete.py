#!/usr/bin/python
import sys
from numpy import *
from scipy import integrate
from matplotlib.pyplot import *

#Define the von mises density of probability
#1st Define Bessel function of order 0
def Io(K):
	f = lambda x : exp(K*cos(x))
	return 1/pi * integrate.quad(f, 0, pi)[0] #[0] means that we just want the result without the error
#2nd Use the Bessel fuunction to define the von mises density of probabilit
#x should be between -pi and +pi
def dp_vonMises(x,mu,K):
	return exp(K*cos(x-mu)) / (2*pi*Io(K))

def baseChange(x, mu):
	return (x-mu)*2*pi

def modified_dp_vonMises(i,mu,K):
	assert 0 <= i and i <= 1
	#res=[]
	#for i in X:
	#	assert 0 <= i and i <= 1
	if 0 <= i and i < mu:
		if mu-0.5 < 0:
			y = dp_vonMises(baseChange(i,mu),0,K) + dp_vonMises(baseChange(i/mu*(mu-0.5), mu), 0,K)*abs(mu-0.5)/mu
		else:
			y = dp_vonMises(baseChange(i, mu), 0,K)
	elif mu <= i and i <= mu + 0.5:
		if mu + 0.5 > 1:
			y = dp_vonMises(baseChange(i,mu),0,K) + dp_vonMises(baseChange((1-i)/(1-mu)*(mu-0.5)+1,mu),0,K)*(mu-0.5)/(1-mu)
		else:
			y = dp_vonMises(baseChange(i,mu),0,K)
	else:
		y = 0
	#	res.append(y*2*pi)
	return y*2*pi
	#return res

#Draw the modified Von Mises density of probability
mu = float(sys.argv[1])
K = float(sys.argv[2])


X1 = linspace(-pi+mu, +pi+mu, 100) # make a list
#Figure with the original von Mises
fig1=figure()
Y1dpVM= dp_vonMises(X1,mu,K)
ax1=fig1.add_subplot(111)
line1=ax1.plot(X1, Y1dpVM, linestyle='-', color=(0.25,0.25,0.25), label="vonMise distribution \n mu=%s and kappa=%s" % (mu, K), lw=2)
handle1, label1 = ax1.get_legend_handles_labels()
ylabel("von Mises distribution")
ax1.fill_between(X1, Y1dpVM, 0, color=(0.75,0.75,0.75))
Y1IntdpVM=[integrate.quad(dp_vonMises, -pi+mu, i, args=(mu,K))[0] for i in X1]
ax2=fig1.add_subplot(111, sharex=ax1, frameon=False)
#ax2=ax1.twinx()
line2=ax2.plot(X1, Y1IntdpVM, linestyle='-.', color='k', label="cumulated probability", lw=2 )
handle2, label2 = ax2.get_legend_handles_labels()
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ylabel("cumulated probability")
#legend((line1, line2), ('von Mises distribution', 'integrated probability'))
#p = Rectangle((0, 0), 1, 1, fc="r")
#legend([p], ["Red Rectangle"])
legend(handle1 + handle2, label1 + label2, loc='upper left')
savefig("1.svg",format='svg')


#Figure with the modified vonMises
X2 = linspace(0, 1, 100) # make a list
#The following allows the function handle lists:
#vmodified_dp_vonMises = vectorize(modified_dp_vonMises)
fig2=figure()
ax1=fig2.add_subplot(111)
vmodified_dp_vonMises = vectorize(modified_dp_vonMises)
Y2dpMVM = vmodified_dp_vonMises(X2,mu,K)
#line1 = ax1.plot(X2,Y2dpMVM, linestyle='-', color=(0.25,0.25,0.25), label="modified von Mises distribution", lw=2)
#handle1, label1 = ax1.get_legend_handles_labels()
#ax1.fill_between(X2, Y2dpMVM, 0, color=(0.75,0.75,0.75))
#ylabel("modified von Mises distribution axis")
Y2IntdpMVM=[integrate.quad(modified_dp_vonMises, 0, i, args=(mu,K))[0] for i in X2]
##ax2 = ax1.twinx()
#ax2=fig2.add_subplot(111, sharex=ax1, frameon=False)
#line2 = ax2.plot(X2,Y2IntdpMVM, linestyle='-.', color='k', label="integrated probability", lw=2)
#handle2, label2 = ax2.get_legend_handles_labels()
#ax2.yaxis.tick_right()
#ax2.yaxis.set_label_position("right")
#ylabel("integration axis")
#legend(handle1 + handle2, label1 + label2)
#savefig("2.svg",format='svg')


def axes_broken_y(axes, upper_frac=0.5, break_frac=0.02, ybounds=None, xlabel=None, ylabel=None):
    """
    Replace the current axes with a set of upper and lower axes.

    The new axes will be transparent, with a breakmark drawn between them.  They
    share the x-axis.  Returns (upper_axes, lower_axes).

    If ybounds=[ymin_lower, ymax_lower, ymin_upper, ymax_upper] is defined,
    upper_frac will be ignored, and the y-axis bounds will be fixed with the
    specified values.
    """
    def breakmarks(axes, y_min, y_max, xwidth=0.008):
        x1, y1, x2, y2 = axes.get_position().get_points().flatten().tolist()
        segment_height = (y_max - y_min) / 3.
        xoffsets = [0, +xwidth, -xwidth, 0]
        yvalues  = [y_min + (i * segment_height) for i in range(4)]
        # Get color of y-axis
        for loc, spine in axes.spines.iteritems():
            if loc  == 'left':
                color = spine.get_edgecolor()
        for x_position in [x1, x2]:
            line = matplotlib.lines.Line2D(
                [x_position + offset for offset in xoffsets], yvalues,
                transform=gcf().transFigure, clip_on=False,
                color=color)
            axes.add_line(line)
    # Readjust upper_frac if ybounds are defined
    if ybounds:
        if len(ybounds) != 4:
            print >> sys.stderr,  "len(ybounds) != 4; aborting..."
            return
        ymin1, ymax1, ymin2, ymax2 = [float(value) for value in ybounds]
        data_height1, data_height2 = (ymax1 - ymin1), (ymax2 - ymin2)
        #upper_frac = data_height2 / (data_height1 + data_height2)
    x1, y1, x2, y2 = axes.get_position().get_points().flatten().tolist()
    width = x2 - x1
    lower_height = (y2 - y1) * ((1 - upper_frac) - 0.5 * break_frac)
    upper_height = (y2 - y1) * (upper_frac - 0.5 * break_frac)
    upper_bottom = (y2 - y1) - upper_height + y1
    lower_axes = matplotlib.pyplot.axes([x1, y1, width, lower_height], axisbg='None')
    upper_axes = matplotlib.pyplot.axes([x1, upper_bottom, width, upper_height],
                          axisbg='None', sharex=lower_axes)
    # Erase the edges between the axes
    for loc, spine in upper_axes.spines.iteritems():
        if loc == 'bottom':
            spine.set_color('none')
    for loc, spine in lower_axes.spines.iteritems():
        if loc == 'top':
            spine.set_color('none')
    upper_axes.get_xaxis().set_ticks_position('top')
    lower_axes.get_xaxis().set_ticks_position('bottom')
    setp(upper_axes.get_xticklabels(), visible=False)
    breakmarks(upper_axes, y1 + lower_height, upper_bottom)
    # Set ylims if ybounds are defined
    if ybounds:
        lower_axes.set_ylim(ymin1, ymax1)
        upper_axes.set_ylim(ymin2, ymax2)
        lower_axes.set_autoscaley_on(False)
        upper_axes.set_autoscaley_on(False)
        upper_axes.yaxis.get_label().set_position((0, 1 - (0.5 / (upper_frac/(1+break_frac)))))
        lower_axes.yaxis.get_label().set_position((0, 0.5 / ((1 - upper_frac)/(1+break_frac))))
    # Make original axes invisible
    axes.set_xticks([])
    axes.set_yticks([])
    #print >> sys.stderr, upper_axes.yaxis.get_label().get_position()
    #print >> sys.stderr, lower_axes.yaxis.get_label().get_position()
    #print >> sys.stderr, axes.yaxis.get_label().get_position()
    #print >> sys.stderr, axes.yaxis.labelpad
    for loc, spine in axes.spines.iteritems():
        spine.set_color('none')
    return upper_axes, lower_axes

def prepare_efficiency(axes, lower_bound=0.69):
    """
    Set up an efficiency figure with breakmarks to indicate a suppressed zero.

    The y-axis limits are set to (lower_bound, 1.0), as appropriate for an
    efficiency plot, and autoscaling is turned off.
    """
    upper_axes, lower_axes = axes_broken_y(axes, upper_frac=0.97)
    lower_axes.set_yticks([])
    upper_axes.set_ylim(lower_bound, 1.)
    upper_axes.set_autoscaley_on(False)
    return upper_axes, lower_axes

print >> sys.stderr, Y2dpMVM[0]

fig3=figure()
ax1=fig3.add_subplot(111)
#ax = plt.axes()
xlim(0,0.5)
upper1, lower1 = axes_broken_y(ax1, ybounds=[0., 10, 150, 250], upper_frac=0.3, break_frac=0.06)
upper1.plot(X2,Y2dpMVM, linestyle='-', color=(0.25,0.25,0.25), label="modified von Mises distribution", lw=2)
xlim(0,0.5)
upper1.fill_between(X2, Y2dpMVM, 0, color=(0.75,0.75,0.75))
xlim(0,0.5)
lower1.plot(X2,Y2dpMVM, linestyle='-', color=(0.25,0.25,0.25), label="modified von Mises distribution", lw=2)
xlim(0,0.5)
lower1.fill_between(X2, Y2dpMVM, 0, color=(0.75,0.75,0.75))
xlim(0,0.5)
handle1, label1 = lower1.get_legend_handles_labels()
ylabel("distribution axis")
xlim(0,0.3)
ax1.spines['bottom'].set_visible(False)

ax2=fig3.add_subplot(111, sharex=ax1, frameon=False)
#upper2, lower2 = axes_broken_y(ax2, ybounds=[0., 10, 150, 250], upper_frac=0.3, break_frac=0.06)
#upper2.plot(X2,Y2IntdpMVM, linestyle='-.', color='k', label="integrated probability", lw=2)
#lower2.plot(X2,Y2IntdpMVM, linestyle='-.', color='k', label="integrated probability", lw=2)
ax2.plot(X2,Y2IntdpMVM, linestyle='-.', color='k', label="cumulated probability", lw=2)
xlim(0,0.5)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ylabel("cumulated probability axis")
xlabel("proportion of the length of the chromosome \n for the slice to be rearranged")
handle2, label2 = ax2.get_legend_handles_labels()
legend(handle1 + handle2, label1 + label2)
xlim(0,0.5)
savefig("3.svg",format='svg')
