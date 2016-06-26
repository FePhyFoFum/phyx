import Tkinter,os,sys
from Tkinter import *


class GUI:
	def __init__(self, master):
		
		self.master = master
		master.title("PHYX GUI")
		
		self.label = Label(master, text="Path to -a").grid(row=0,column=0)
		#self.label.pack()
		
		global a
		a = Tkinter.StringVar()
		self.a_seq = Entry(master, bd=5, textvariable=a).grid(row=0,column=1)
		#self.aa_seq.pack(side = LEFT)

		self.label = Label(master, text="Path to -b").grid(row=0,column=2)
		#self.label.pack(side = LEFT)

		global b
		b = Tkinter.StringVar()
		self.b_seq = Entry(master, bd=5, textvariable=b).grid(row=0,column=3)
		#self.b_seq.pack(side = LEFT)
		
		self.label = Label(master, text="Path to -c").grid(row=0,column=4)
		#self.label.pack(side = LEFT)
		
		global c
		c = Tkinter.StringVar()
		self.c_seq = Entry(master, bd=5, textvariable=c).grid(row=0,column=5)
		#self.c_seq.pack(side = LEFT)
		
		self.label = Label(master, text="Path to -d").grid(row=0,column=6)
		#self.label.pack(side = LEFT)
		
		global d
		d = Tkinter.StringVar()
		self.d_seq = Entry(master, bd=5, textvariable=d).grid(row=0,column=7)
		#self.d_seq.pack(side = LEFT)
		
		self.label = Label(master, text="Path to -e").grid(row=1,column=0)		
		global e
		e = Tkinter.StringVar()
		self.e_seq = Entry(master, bd=5, textvariable=e).grid(row=1,column=1)
		
		self.label = Label(master, text="Path to -f").grid(row=1,column=2)		
		global f
		f = Tkinter.StringVar()
		self.f_seq = Entry(master, bd=5, textvariable=f).grid(row=1,column=3)
		
		self.label = Label(master, text="Path to -g").grid(row=1,column=4)		
		global g
		g = Tkinter.StringVar()
		self.g_seq = Entry(master, bd=5, textvariable=g).grid(row=1,column=5)
		
		self.label = Label(master, text="Path to -h").grid(row=1,column=6)		
		global h
		h = Tkinter.StringVar()
		self.h_seq = Entry(master, bd=5, textvariable=h).grid(row=1,column=7)
		
		self.label = Label(master, text="Path to -i").grid(row=2,column=0)		
		global i
		i = Tkinter.StringVar()
		self.i_seq = Entry(master, bd=5, textvariable=i).grid(row=2,column=1)
		
		self.label = Label(master, text="Path to -j").grid(row=2,column=2)		
		global j
		j = Tkinter.StringVar()
		self.j_seq = Entry(master, bd=5, textvariable=j).grid(row=2,column=3)
		
		self.label = Label(master, text="Path to -k").grid(row=2,column=4)		
		global k
		k = Tkinter.StringVar()
		self.k_seq = Entry(master, bd=5, textvariable=k).grid(row=2,column=5)
		
		self.label = Label(master, text="Path to -l").grid(row=2,column=6)		
		global l
		l = Tkinter.StringVar()
		self.l_seq = Entry(master, bd=5, textvariable=l).grid(row=2,column=7)
		
		self.label = Label(master, text="Path to -m").grid(row=3,column=0)		
		global m
		m = Tkinter.StringVar()
		self.m_seq = Entry(master, bd=5, textvariable=m).grid(row=3,column=1)
		
		self.label = Label(master, text="Path to -n").grid(row=3,column=2)		
		global n
		n = Tkinter.StringVar()
		self.n_seq = Entry(master, bd=5, textvariable=n).grid(row=3,column=3)
		
		self.label = Label(master, text="Path to -o").grid(row=3,column=4)		
		global o
		o = Tkinter.StringVar()
		self.o_seq = Entry(master, bd=5, textvariable=o).grid(row=3,column=5)
		
		self.label = Label(master, text="Path to -p").grid(row=3,column=6)		
		global p
		p = Tkinter.StringVar()
		self.p_seq = Entry(master, bd=5, textvariable=p).grid(row=3,column=7)
		
		self.label = Label(master, text="Path to -q").grid(row=4,column=0)		
		global q
		q = Tkinter.StringVar()
		self.q_seq = Entry(master, bd=5, textvariable=q).grid(row=4,column=1)
		
		self.label = Label(master, text="Path to -r").grid(row=4,column=2)		
		global r
		r = Tkinter.StringVar()
		self.r_seq = Entry(master, bd=5, textvariable=r).grid(row=4,column=3)
		
		self.label = Label(master, text="Path to -s").grid(row=4,column=4)		
		global s
		s = Tkinter.StringVar()
		self.s_seq = Entry(master, bd=5, textvariable=s).grid(row=4,column=5)
		
		self.label = Label(master, text="Path to -t").grid(row=4,column=6)		
		global t
		t = Tkinter.StringVar()
		self.t_seq = Entry(master, bd=5, textvariable=t).grid(row=4,column=7)
		
		self.label = Label(master, text="Path to -u").grid(row=5,column=0)		
		global u
		u = Tkinter.StringVar()
		self.u_seq = Entry(master, bd=5, textvariable=u).grid(row=5,column=1)
		
		self.label = Label(master, text="Path to -v").grid(row=5,column=2)		
		global v
		v = Tkinter.StringVar()
		self.v_seq = Entry(master, bd=5, textvariable=v).grid(row=5,column=3)
		
		self.label = Label(master, text="Path to -w").grid(row=5,column=4)		
		global w
		w = Tkinter.StringVar()
		self.w_seq = Entry(master, bd=5, textvariable=w).grid(row=5,column=5)
		
		self.label = Label(master, text="Path to -x").grid(row=5,column=6)		
		global x
		x = Tkinter.StringVar()
		self.x_seq = Entry(master, bd=5, textvariable=x).grid(row=5,column=7)
		
		self.label = Label(master, text="Path to -y").grid(row=6,column=0)		
		global y
		y = Tkinter.StringVar()
		self.y_seq = Entry(master, bd=5, textvariable=y).grid(row=6,column=1)
		
		self.label = Label(master, text="Path to -z").grid(row=6,column=2)		
		global z
		z = Tkinter.StringVar()
		self.z_seq = Entry(master, bd=5, textvariable=z).grid(row=6,column=3)
		"""
		self.label = Label(master, text="Path to -n")
		self.label.pack(padx=10)
		
		global test2
		test2 = Tkinter.StringVar()
		self.nuc_seq = Entry(master, bd=5, textvariable=test2)
		self.nuc_seq.pack(padx=10)
		
		self.label = Label(master, text="Outfile -o")
		self.label.pack()
		
		global test3
		test3 = ""
		test3 = Tkinter.StringVar()
		self.out_seq = Entry(master, bd=5, textvariable=test3)
		self.out_seq.pack()
		
		self.label = Label(master, text="help -h")
		self.label.pack()
		
		global help_me
		help_me = ""
		help_me = Tkinter.StringVar()
		self.help_seq = Entry(master, bd=5, textvariable=help_me)
		self.help_seq.pack()
		"""
		self.aa_button = Button(master, text="Run pxaatocdn", command=self.aa).grid(row=7,column=0)
		#self.aa_button.pack()
		
		self.bdfit_button = Button(master, text="Run pxbdfit", command=self.bdfit).grid(row=7,column=1)
		
		self.bdsim_button = Button(master, text="Run pxbdsim", command=self.bdsim).grid(row=7,column=2)
		
		self.boot_button = Button(master, text="Run pxboot", command=self.boot).grid(row=7,column=3)
		
		self.bp_button = Button(master, text="Run pxbp", command=self.bp).grid(row=7,column=4)
		
		self.bpseq_button = Button(master, text="Run pxbpseq", command=self.bpseq).grid(row=7,column=5)

		self.cat_button = Button(master, text="Run pxcat", command=self.cat).grid(row=7,column=6)

		self.clsq_button = Button(master, text="Run pxclsq", command=self.clsq).grid(row=7,column=7)

		self.conseq_button = Button(master, text="Run pxconseq", command=self.conseq).grid(row=8,column=0)

		self.contrates_button = Button(master, text="Run pxcontrates", command=self.contrates).grid(row=8,column=1)

		self.fqfilt_button = Button(master, text="Run pxfqfilt", command=self.fqfilt).grid(row=8,column=2)

		self.lsseq_button = Button(master, text="Run pxlsseq", command=self.lsseq).grid(row=8,column=3)

		self.lstr_button = Button(master, text="Run pxlstr", command=self.lstr).grid(row=8,column=4)
		
		self.mrca_button = Button(master, text="Run pxmrca", command=self.mrca).grid(row=8,column=5)

		self.mrcacut_button = Button(master, text="Run pxmrcacut", command=self.mrcacut).grid(row=8,column=6)

		self.mrcaname_button = Button(master, text="Run pxmrcaname", command=self.mrcaname).grid(row=8,column=7)

		self.nj_button = Button(master, text="Run pxnj", command=self.nj).grid(row=9,column=0)

		self.nni_button = Button(master, text="Run pxnni", command=self.nni).grid(row=9,column=1)

		self.nw_button = Button(master, text="Run pxnw", command=self.nw).grid(row=9,column=2)

		self.recode_button = Button(master, text="Run pxrecode", command=self.recode).grid(row=9,column=3)

		self.revcomp_button = Button(master, text="Run pxrevcomp", command=self.revcomp).grid(row=9,column=4)

		self.rms_button = Button(master, text="Run pxrms", command=self.rms).grid(row=9,column=5)

		self.rmt_button = Button(master, text="Run pxrmt", command=self.rmt).grid(row=9,column=6)

		self.rr_button = Button(master, text="Run pxrr", command=self.rr).grid(row=9,column=7)

		self.seqgen_button = Button(master, text="Run pxseqgen", command=self.seqgen).grid(row=10,column=0)

		self.stofa_button = Button(master, text="Run pxstofa", command=self.stofa).grid(row=10,column=1)

		self.stonex_button = Button(master, text="Run pxstonex", command=self.stonex).grid(row=10,column=2)

		self.stophy_button = Button(master, text="Run pxstophy", command=self.stophy).grid(row=10,column=3)

		self.strec_button = Button(master, text="Run pxstrec", command=self.strec).grid(row=10,column=4)

		self.sw_button = Button(master, text="Run pxsw", command=self.sw).grid(row=10,column=5)

		self.tlate_button = Button(master, text="Run pxtlate", command=self.tlate).grid(row=10,column=6)

		self.upgma_button = Button(master, text="Run pxupgma", command=self.upgma).grid(row=10,column=7)

		#self.seqgen_button.pack()
		
		self.h_button = Button(master, text="print help menu", command=self.h).grid(row=11,column=3)
		#self.h_button.pack()
		#global help_on
		#help_on = ""
		#help_on = Tkinter.StringVar()
		
		self.close_button = Button(master, text="Close program", command=master.quit).grid(row=11,column=4)
		#self.close_button.pack()
	def aa(self):
		if (h.get() != ""):
			cmd = "~/pxaatocdn -h"
		else: 
			cmd = "~/pxaatocdn -a " + a.get() + " -n " + n.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def bdfit(self):
		if (h.get() != ""):
			cmd = "~/pxbdfit -h"
		else: 
			cmd = "~/pxbdfit -t " + t.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def bdsim(self):
		if (h.get() != ""):
			cmd = "~/pxbdsim -h"
		else: 
			cmd = "~/pxbdsim "
			if (e.get() != ""):
				cmd += " -e " + e.get()
			if (t.get() != ""):
				cmd += " -t " + t.get()
			if (b.get() != ""):
				cmd += " -b " + b.get()
			if (d.get() != ""):
				cmd += " -d " + d.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (s.get() != ""):
				cmd += " -s " + s.get()
			if (x.get() != ""):
				cmd += " -x " + x.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def boot(self):
		if (h.get() != ""):
			cmd = "~/pxboot -h"
		else: 
			cmd = "~/pxboot "
			if (s.get() != ""):
				cmd += " -s " + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
			if (p.get() != ""):
				cmd += " -p " + p.get()
			if (f.get() != ""):
				cmd += " -f " + f.get()
			if (x.get() != ""):
				cmd += " -x " + x.get()
		os.system(cmd)
	def bp(self):
		if (h.get() != ""):
			cmd = "~/pxbp -h"
		else: 
			cmd = "~/pxbp -t" + t.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def bpseq(self):
		if (h.get() != ""):
			cmd = "~/pxbpseq -h"
		else: 
			cmd = "~/pxbpseq -t" + t.get() + " -s " + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def cat(self):
		if (h.get() != ""):
			cmd = "~/pxcat -h"
		else: 
			cmd = "~/pxcat -s" + s.get()
			if (p.get() != ""):
				cmd += " -p " + p.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def clsq(self):
		if (h.get() != ""):
			cmd = "~/pxclsq -h"
		else: 
			cmd = "~/pxclsq -s" + s.get()
			if (p.get() != ""):
				cmd += " -p " + p.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def conseq(self):
		if (h.get() != ""):
			cmd = "~/pxconseq -h"
		else: 
			cmd = "~/pxconseq -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def contrates(self):
		if (h.get() != ""):
			cmd = "~/pxcontrates -h"
		else: 
			cmd = "~/pxcontrates -s" + s.get()
			if (c.get() != ""):
				cmd += " -c " + c.get()
			if (t.get() != ""):
				cmd += " -t " + t.get()
			if (a.get() != ""):
				cmd += " -a " + a.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def fqfilt(self):
		if (h.get() != ""):
			cmd = "~/pxfqfilt -h"
		else: 
			cmd = "~/pxfqfilt -s" + s.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (s.get() != ""):
				cmd += " -s " + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def lsseq(self):
		if (h.get() != ""):
			cmd = "~/pxlsseq -h"
		else: 
			cmd = "~/pxlsseq -s" + s.get()
			if (a.get() != ""):
				cmd += " -a "
			if (p.get() != ""):
				cmd += " -p "
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def lstr(self):
		if (h.get() != ""):
			cmd = "~/pxlstr -h"
		else: 
			cmd = "~/pxlstr -t" + t.get()
			if (a.get() != ""):
				cmd += " -a "
			if (r.get() != ""):
				cmd += " -r "
			if (b.get() != ""):
				cmd += " -b "
			if (u.get() != ""):
				cmd += " -u "
			if (n.get() != ""):
				cmd += " -n "
			if (l.get() != ""):
				cmd += " -l "
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def mrca(self):
		if (h.get() != ""):
			cmd = "~/pxmrca -h"
		else: 
			cmd = "~/pxmrca -t" + t.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def mrcacut(self):
		if (h.get() != ""):
			cmd = "~/pxmrcacut -h"
		else: 
			cmd = "~/pxmrcacut -t" + t.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def mrcaname(self):
		if (h.get() != ""):
			cmd = "~/pxmrcaname -h"
		else: 
			cmd = "~/pxmrcaname -t" + t.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def nj(self):
		if (h.get() != ""):
			cmd = "~/pxnj -h"
		else: 
			cmd = "~/pxnj -s" + s.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def nni(self):
		if (h.get() != ""):
			cmd = "~/pxnni -h"
		else: 
			cmd = "~/pxnni -t" + t.get()
			if (x.get() != ""):
				cmd += " -x " + x.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def nw(self):
		if (h.get() != ""):
			cmd = "~/pxnw -h"
		else: 
			cmd = "~/pxnw -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
			if (a.get() != ""):
				cmd += " -a " + a.get()
			if (t.get() != ""):
				cmd += " -t " + t.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (v.get() != ""):
				cmd += " -v " + v.get()
		os.system(cmd)
	def recode(self):
		if (h.get() != ""):
			cmd = "~/pxrecode -h"
		else: 
			cmd = "~/pxrecode -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def revcomp(self):
		if (h.get() != ""):
			cmd = "~/pxrevcomp -h"
		else: 
			cmd = "~/pxrevcomp -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def rms(self):
		if (h.get() != ""):
			cmd = "~/pxrms -h"
		else: 
			cmd = "~/pxrms -s" + s.get()
			if (r.get() != ""):
				cmd += " -r " + r.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def rmt(self):
		if (h.get() != ""):
			cmd = "~/pxrmt -h"
		else: 
			cmd = "~/pxrmt -t" + t.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (f.get() != ""):
				cmd += " -f " + f.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def rr(self):
		if (h.get() != ""):
			cmd = "~/pxrr -h"
		else: 
			cmd = "~/pxrr -t" + t.get()
			if (f.get() != ""):
				cmd += " -f " + f.get()
			if (u.get() != ""):
				cmd += " -u " + u.get()
			if (s.get() != ""):
				cmd += " -s " + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def seqgen(self):
		if (h.get() != ""):
			cmd = "~/pxseqgen -h"
		else: 
			cmd = "~/pxseqgen -t" + t.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
			if (l.get() != ""):
				cmd += " -l " + l.get()
			if (b.get() != ""):
				cmd += " -b " + b.get()
			if (g.get() != ""):
				cmd += " -g " + g.get()
			if (r.get() != ""):
				cmd += " -r " + r.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (x.get() != ""):
				cmd += " -x " + x.get()
			if (a.get() != ""):
				cmd += " -a "
			if (p.get() != ""):
				cmd += " -p "
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (k.get() != ""):
				cmd += " -k " + k.get()
		os.system(cmd)
	def stofa(self):
		if (h.get() != ""):
			cmd = "~/pxstofa -h"
		else: 
			cmd = "~/pxstofa -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def stonex(self):
		if (h.get() != ""):
			cmd = "~/pxstonex -h"
		else: 
			cmd = "~/pxstonex -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def stophy(self):
		if (h.get() != ""):
			cmd = "~/pxstophy -h"
		else: 
			cmd = "~/pxstophy -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def strec(self):
		if (h.get() != ""):
			cmd = "~/pxstrec -h"
		else: 
			cmd = "~/pxstrec -d" + d.get()
			if (w.get() != ""):
				cmd += " -w "
			if (z.get() != ""):
				cmd += " -z "
			if (t.get() != ""):
				cmd += " -t " + t.get()
			if (c.get() != ""):
				cmd += " -c " + c.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (a.get() != ""):
				cmd += " -a " + a.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (p.get() != ""):
				cmd += " -p " + p.get()
			if (l.get() != ""):
				cmd += " -l " + l.get()
		os.system(cmd)
	def sw(self):
		if (h.get() != ""):
			cmd = "~/pxsw -h"
		else: 
			cmd = "~/pxsw -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
			if (a.get() != ""):
				cmd += " -a " + a.get()
			if (t.get() != ""):
				cmd += " -t " + t.get()
			if (m.get() != ""):
				cmd += " -m " + m.get()
			if (n.get() != ""):
				cmd += " -n " + n.get()
			if (v.get() != ""):
				cmd += " -v " + v.get()
		os.system(cmd)
	def tlate(self):
		if (h.get() != ""):
			cmd = "~/pxtlate -h"
		else: 
			cmd = "~/pxtlate -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def upgma(self):
		if (h.get() != ""):
			cmd = "~/pxupgma -h"
		else: 
			cmd = "~/pxtlate -s" + s.get()
			if (o.get() != ""):
				cmd += " -o " + o.get()
		os.system(cmd)
	def h(self):
		print "========================================================="
		print "This is GUI for PHYX designed to make input easier"
		print "Type the location of the file or value into the box"
		print "For help on a program, type anything into the help box"
		print "and click the program you want help with. Not all options are"
		print "needed for all programs so only fill in boxes used for "
		print "the program of interest. Anything extra filled in will be" 
		print "ignored if it is not an option for the program"
		print "If the programs option does not require an argument"
		print "simply type something into the box and it will activate"
		print "that argument"
		print "Copyright (C) 2015 FePhyFoFum\nLicense GPLv2\nwritten by Joseph F. Walker, Joseph W. Brown, Stephen A. Smith (blackrim)"
		print "========================================================="


root = Tk()
my_gui = GUI(root)
root.mainloop()
