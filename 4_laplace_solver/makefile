ALL=4_1_laplace_solver 4_2_1_laplace_solver 4_2_2_laplace_solver 4_3_laplace_solver 4_4_laplace_solver

THREADS=24

BERKELEY_UPC=0
CRAY_UPC=1

NETWORK=smp
BERKELEY_FLAGS=-pthreads

UPCRUN=upcrun

ifeq ($(strip $(BERKELEY_UPC)),1)
    UPCFLAGS := -network=$(NETWORK) -pthreads -T=$(THREADS) $(BERKELEY_FLAGS)
    UPCC = upcc
else
ifeq ($(strip $(CRAY_UPC)),1)
	UPCFLAGS := -h upc -X $(THREADS)
	UPCC = cc
else
    UPCFLAGS := -Wall -T=$(THREADS)
    UPCC = upc
endif
endif

all: $(ALL)

%.o: %.upc
	+ $(UPCC) -c $(UPCFLAGS) $< -o $@

%: %.o
	+ $(UPCC) $(UPCFLAGS) $< -o $@

run%: %
	+ $(UPCRUN) -n $(THREADS) $<

clean:
	$(RM) -rf $(ALL) *_pthread-link job_upc_ex*.e*  job_upc_ex*.o*