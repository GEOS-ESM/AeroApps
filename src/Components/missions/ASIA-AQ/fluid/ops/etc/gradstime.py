import datetime as dt

class GradsTime(object):

    def __init__(self, idate, itime, shift_dt=0):

        idate = str(idate)
        itime = '%06d'%(itime,)

        self.idt = dt.datetime.strptime(idate+itime,'%Y%m%d%H%M%S')
        self.idt += dt.timedelta(hours=shift_dt)

    def strftime(self, s, tau):

        return self.strvtime(self.stritime(s), tau)

    def strvtime(self, s, tau):

        vdt = self.idt + dt.timedelta(hours=tau)

        tokens = { '%y2' : vdt.strftime('%y'),
                   '%y4' : vdt.strftime('%Y'),
                   '%m1' : str(vdt.month),
                   '%m2' : vdt.strftime('%m'),
                   '%mc' : vdt.strftime('%b'),
                   '%d1' : str(vdt.day),
                   '%d2' : vdt.strftime('%d'),
                   '%d2' : vdt.strftime('%d'),
                   '%h1' : str(vdt.hour),
                   '%h2' : vdt.strftime('%H'),
                   '%h3' : '%03d'%(tau,),
                   '%f2' : '%02d'%(tau,),
                   '%f3' : '%03d'%(tau,),
                   '%n2' : vdt.strftime('%M'),
                   '%j3' : vdt.strftime('%j')
                 }

        for k,v in tokens.iteritems(): s = s.replace(k,v)

        return s

    def stritime(self, s, tau=None):

        tokens = { '%iy2' : self.idt.strftime('%y'),
                   '%iy4' : self.idt.strftime('%Y'),
                   '%im1' : str(self.idt.month),
                   '%im2' : self.idt.strftime('%m'),
                   '%imc' : self.idt.strftime('%b'),
                   '%id1' : str(self.idt.day),
                   '%id2' : self.idt.strftime('%d'),
                   '%id2' : self.idt.strftime('%d'),
                   '%ih1' : str(self.idt.hour),
                   '%ih2' : self.idt.strftime('%H'),
                   '%ih3' : '%03d'%(self.idt.hour,),
                   '%in2' : self.idt.strftime('%M'),
                   '%ij3' : self.idt.strftime('%j')
                 }

        for k,v in tokens.iteritems(): s = s.replace(k,v)

        return s
