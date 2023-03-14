include make.inc

LIBOBJS := $(OBJDIR)/ab2rht.o $(OBJDIR)/utils.o $(OBJDIR)/rht2ht.o $(OBJDIR)/blockreflectors.o

all: $(BINDIR)/profilereduction

# 
# Lib objs
# 

$(OBJDIR)/utils.o: $(SRCDIR)/utils.f90 make.inc
	$(FC) $(FFLAGS) -o $@ -J$(INCDIR) -I$(INCDIR) -c $<

$(OBJDIR)/dgghd3.o: $(SRCDIR)/dgghd3.f $(OBJDIR)/utils.o make.inc
	$(FC) $(FFLAGS) -o $@ -J$(INCDIR) -I$(INCDIR) -c $<

$(OBJDIR)/dgghd4.o: $(SRCDIR)/dgghd4.f $(OBJDIR)/utils.o make.inc
	$(FC) $(FFLAGS) -o $@ -J$(INCDIR) -I$(INCDIR) -c $<

$(OBJDIR)/ab2rht.o: $(SRCDIR)/ab2rht.f90 $(OBJDIR)/utils.o $(OBJDIR)/blockreflectors.o make.inc
	$(FC) $(FFLAGS) -o $@ -J$(INCDIR) -I$(INCDIR) -c $<

$(OBJDIR)/rht2ht.o: $(SRCDIR)/rht2ht.f90 $(OBJDIR)/utils.o make.inc
	$(FC) $(FFLAGS) -o $@ -J$(INCDIR) -I$(INCDIR) -c $<

$(OBJDIR)/blockreflectors.o: $(SRCDIR)/blockreflectors.f90 $(OBJDIR)/utils.o make.inc
	$(FC) $(FFLAGS) -o $@ -J$(INCDIR) -I$(INCDIR) -c $<


# 
# Executables
# 

$(BINDIR)/profilereduction: $(SRCDIR)/profilereduction.f90 $(LIBOBJS) make.inc
	$(FC) $(FFLAGS) -o $@ $<  $(LIBOBJS) -I$(INCDIR) $(LIBS) -lstdc++

clean:
	rm -f ./*/*.mod $(BINDIR)/* $(INCDIR)/* $(OBJDIR)/*