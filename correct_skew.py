#!/usr/bin/env python

from pylab import *
from sympy import *
from sys import *
from numpy import *
from re import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats.mstats import mquantiles

def line(slope_,y_intercept):
	"Returns a line object from a slope and y_intercept."
	return Line(Point(0,y_intercept),slope=slope_)

def xy_coord(l1,l2):
	"Returns the x,y-coordinates of intersection between l1 and l2"
	p=l1.intersection(l2)
	return p[0]

def analyze_delays(delays,delays_no_skew, src_dest,arrivals,alpha):
	"""
	Plot delays before and after skew corrections
	"""
	diff = [ d-d_ for d in delays_no_skew for d_ in delays ]

	pp = PdfPages("Clock_skew_correction.pdf")

	fig1 = figure()
	ax = fig1.add_subplot(111)
	hist(diff,50,facecolor='green')
	ax.set_title("Histogram of corrections"+"; src:dest:size"+ src_dest+"\nclock ratio:"+repr(alpha))
	ax.set_xlabel("corrections (ns)")
	ax.set_ylabel("count")
	pp.savefig(fig1)

	fig2 = figure()
	ax1 = fig2.add_subplot(111)
	ax1.scatter(arrivals,delays,c='b',label='with-skew')
	ax1.scatter(arrivals,delays_no_skew,c='r',label='no-skew')
	legend(loc='upper left')
	pp.savefig(fig2)

	pp.close()


def preprocess_prv(prv,method=False,by_node=False):
	"""
	1) Parses the paraver files
	2) Computes delays using method1 (precv-lrecv), or method2 (precv-lsend)
	3) Group comms by nodes (if by_node is true ) or by comms (if by_node is false)
	4) Invokes estimate_ratio
	comm_fields:map

	3:object_send:lsend:psend:object_recv:lrecv:precv:size:tag
	object_send: [cpu_send_id:ptask_send_id:task_send_id:thread_send_id]
	object_recv: [cpu_recv_id:ptask_recv_id:task_recv_id:thread_recv_id]

	dep_t_s1=0
	dep_t_si: time duration between the first and the ith packets' departures
				consistent with the sender's clock (Cs).

	arriv_t_r1=0
	arriv_t_ri: time duration between the first and the ith packets' arrivals
				at the receiver consistent with the receiver's clock (Cr).
							key,value
	For each packet {from_rank:to_rank:size,[delay_f,dep_t_si,arriv_t_ri]}
	"""
	comms_by_rank_size={}

	with open(prv_file,'r') as PRV:
		for line in PRV:
			line=line.strip('\n')
			comm=compile("^3:.+")
			init=[]
			if comm.match(line):
				comm_fields=line.split(':')
				key=comm_fields[1]+':'+comm_fields[7]+':'+comm_fields[-2]
				delay_f=int(comm_fields[12])-int(comm_fields[5])

				if key in comms_by_rank_size.keys():
					dep_t_si = int(comm_fields[5]) - int(comms_by_rank_size[key][0][-2])
					arriv_t_ri = int(comm_fields[12]) - int(comms_by_rank_size[key][0][-1])
					arriv_tr = int(comm_fields[12])
					comms_by_rank_size[key].append((delay_f,dep_t_si,arriv_tr,arriv_t_ri,arriv_t_ri-dep_t_si))
				else:
					init.append((delay_f,0,int(comm_fields[12]),0,0,int(comm_fields[5]),int(comm_fields[12])))
					comms_by_rank_size[key]=init

	pp = PdfPages("Clock_skew_correction.pdf")
	#for key in comms_by_rank_size.iterkeys():
	for key in comms_by_rank_size.keys()[:50]:
		delays=[ rec[4] for rec in comms_by_rank_size[key] ]
		send_s=[ rec[1] for rec in comms_by_rank_size[key] ]
		arrival_r=[ rec[2] for rec in comms_by_rank_size[key] ]
		if len(delays) < 10:
			continue
		alpha,beta = estimate_ratio(delays,send_s,len(delays))

		delays_no_skew = [ delays[i]-((alpha-1)*send_s[i])+beta for i in xrange(0,len(delays))]
		diff = [ delays_no_skew[i] - delays[i] for i in xrange(0,len(delays)) ]

		print "len(delays_no_skew):" + repr(len(delays_no_skew))
		print "len(arrival_r):" + repr(len(arrival_r))
		#analyze_delays(delays,delays_no_skew,key,arrivals,alpha)
		#print "alpha: " + repr(alpha) + ", beta: " + repr(beta)

		Qu_1,median,Qu_3 = mquantiles(diff)
		min_,mean_,max_ = min(diff),mean(diff),max(diff)


		fig1 = figure()
		ax = fig1.add_subplot(111)
		hist(diff,50,color='green')
		ax.set_title("Distribution of corrections"+
					"; src:dest:size :: "+
					key +
					"\nclock ratio:"+repr(float64(alpha))+ "total count:"+repr(len(delays)))
					#"\n; min,1_Qu,median,mean,3_Qu,max"+repr(min_)+","+repr(Qu_1)+","+repr(median)+"\n,"+repr(mean_)+","+repr(Qu_3)+","+repr(max_))
		ax.set_xlabel("corrections (ns)")
		ax.set_ylabel("count")
		pp.savefig(fig1)

		fig2 = figure()
		ax1 = fig2.add_subplot(111)
		ax1.scatter(float64(arrival_r),float64(delays),c='b',label='with-skew')
		ax1.scatter(float64(arrival_r),float64(delays_no_skew),c='r',label='no-skew')
		ax1.set_title("Delays: before and after skew correction:\n src:dest:size :: "+key+"clock ratio: "+repr(float64(alpha)))
		ax1.set_xlabel("Arrival time")
		ax1.set_ylabel("delays(ns)")
		legend(loc='upper left')
		pp.savefig(fig2)

	pp.close()

def estimate_ratio(delays,send_s,N):
	alpha,beta=0,0
	k,n=2,array([0,1,2])
	j=0
	d,send_t=asarray(delays),asarray(send_s)
	for i in xrange(3,N):
		j=k
		while j > 2:
			X= xy_coord( line( send_t[i], -1*d[i] ), line( send_t[n[j]], -1*d[n[j]] ) ).x
			X_= xy_coord( line(send_t[n[j]], -1*d[n[j]]), line( send_t[n[j-1]], -1*d[n[j-1]]) ).x
			if X > X_:
				break
			j=j-1
		k=j+1
		n=insert(n,k,i)
	opts=sum(send_t)/N
	for i in xrange(1,k-1):
		if send_t[n[i]] < opts and opts < send_t[n[i+1]]:
			alpha = xy_coord( line(send_t[n[i]],-1*d[n[i]]), line(send_t[n[i+1]], -1*d[n[i+1]]) ).x + 1
			beta = xy_coord( line(send_t[n[i]],-1*d[n[i]]), line(send_t[n[i+1]], -1*d[n[i+1]]) ).y
			break
	return alpha,beta

if __name__ == '__main__':
	prv_file=sys.argv[1]
	preprocess_prv(prv_file)
