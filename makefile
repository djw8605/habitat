LINK_TARGET = habitat.exe


OBJS =  \
 ecologst.o \
 glider.o \
 habitat.o \
 landscpe.o \
 linked.o \
 mt19937.o \
 queue.o \
 util98.o 


REBUILDABLES = $(OBJS) $(LINK_TARGET)

clean : 	
	rm -f $(REBUILDABLES)	
	echo Clean done
 
all : $(LINK_TARGET)	
	echo All done 

$(LINK_TARGET) : $(OBJS)	
	g++ -fopenmp -g -o $@ $^ 

%.o : %.cpp	
	g++ -fopenmp -g -o $@ -c $< 

ecologst.o : ecologst.h glider.h habitat.h landscpe.h linked.h queue.h util98.h
glider.o : glider.h habitat.h landscpe.h linked.h queue.h util98.h
habitat.o : ecologst.h glider.h habitat.h landscpe.h linked.h queue.h util98.h
landscpe.o : habitat.h landscpe.h util98.h
linked.o : habitat.h linked.h util98.h
mt19937.o : mt19937.h 
queue.o : habitat.h queue.h util98.h
util98.o : util98.h mt19937.h


 
