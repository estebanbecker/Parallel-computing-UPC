ALL=Conversion_table

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