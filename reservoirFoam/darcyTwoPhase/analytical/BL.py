#!/usr/bin/env python

class classRockFluid:
    Swi=0.
    Sor=0.
    nw=1.
    no=1.
    krwe=1.0
    krne=1.0
    muw=1e-3
    mun=1e-3

    def krw(self,Sw):
        return self.krwe*((Sw-self.Swi)/(1.-self.Swi-self.Sor))**nw

    def krn(self,Sw):
        return self.krne*((1-self.Sor-Sw)/(1.-self.Swi-self.Sor))**nw

    def Fw(self,Sw):
        return self.krw(Sw)/self.muw / (self.krw(Sw)/self.muw + self.krn(Sw)/self.mun)
