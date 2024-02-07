# Makefile Projet Modelisation pour les energies

# Compiler used
CC = gfortran-11

# Options de compilation
FFLAGS = -O0

# Chemin ou sont stockes les fichiers binaires
ODIR = ../bin
# Chemin ou sera stocke le fichier executable
EXEDIR = ..

# Liste des fichiers binaires
OBJS = $(ODIR)/main.o $(ODIR)/ondelettes.o $(ODIR)/leastsquares.o $(ODIR)/initialisation.o $(ODIR)/numerics.o
# Liste des fichiers mod
MODS = ondelettes.mod leastsquares.mod initialisation.mod numerics.mod

# Librairies used
LIBS = -llapack -lblas

# nom de l'executable
PROG = pb_inverse

# regle de compilation pour le fichier executable final
$(PROG) : $(OBJS)
	$(CC) $(FFLAGS) $(OBJS) $(LIBS) -o $(EXEDIR)/$(PROG)

# regle clean pour supprimer les fichiers binaires et executable
clean:
	rm $(ODIR)/*.o $(EXEDIR)/$(PROG)

# Regles construction des fichiers binaires
$(ODIR)/main.o: main.f90 ondelettes.f90 leastsquares.f90 initialisation.f90 numerics.f90
	$(CC) -c $(FFLAGS) main.f90 -o $(ODIR)/main.o
$(ODIR)/numerics.o: numerics.f90 
	$(CC) -c $(FFLAGS) numerics.f90 -o $(ODIR)/numerics.o
$(ODIR)/ondelettes.o: ondelettes.f90 numerics.f90
	$(CC) -c $(FFLAGS) ondelettes.f90 -o $(ODIR)/ondelettes.o
$(ODIR)/leastsquares.o: leastsquares.f90 numerics.f90
	$(CC) -c $(FFLAGS) leastsquares.f90 -o $(ODIR)/leastsquares.o
$(ODIR)/initialisation.o: initialisation.f90 initialisation.f90
	$(CC) -c $(FFLAGS) initialisation.f90 -o $(ODIR)/initialisation.o
