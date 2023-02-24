.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean))) #############

# CFLAGS := -Wall -O3 -flto
CFLAGS := -Wall -Og -g
CFLAGS += -fmax-errors=3 -Iinclude
# CFLAGS += -DNDEBUG
CPPSTD := c++17

# generate .d files during compilation
DEPFLAGS = -MT $@ -MMD -MP -MF .build/$*.d

FIND_MAIN := \
  find src -type f -regex '.*\.cc?$$' \
  | xargs grep -l '^\s*int\s\+main\s*(' \
  | sed 's:^src/\(.*\)\.c\+$$:bin/\1:'
EXE := $(shell $(FIND_MAIN))

all: $(EXE)

bin/binner: LDFLAGS += -static -static-libgcc -static-libstdc++

ifneq (, $(shell which root-config))
ROOT_CFLAGS  := $(shell root-config --cflags | sed 's/ -std=c++[^ ]\+ / /')
ROOT_LDFLAGS := $(shell root-config --ldflags)
ROOT_LDLIBS  := $(shell root-config --libs)

.build/make_vars.o: CPPSTD = c++20
.build/make_vars.o: CFLAGS += $(ROOT_CFLAGS)
bin/make_vars: LDFLAGS += $(ROOT_LDFLAGS)
bin/make_vars: LDLIBS  += $(ROOT_LDLIBS)

.build/convert_mxaods.o: CPPSTD = c++20
.build/convert_mxaods.o: CFLAGS += $(ROOT_CFLAGS)
bin/convert_mxaods: LDFLAGS += $(ROOT_LDFLAGS)
bin/convert_mxaods: LDLIBS  += $(ROOT_LDLIBS)

.build/varlist.o: CPPSTD = c++20
.build/varlist.o: CFLAGS += $(ROOT_CFLAGS)
bin/varlist: LDFLAGS += $(ROOT_LDFLAGS)
bin/varlist: LDLIBS  += $(ROOT_LDLIBS)
endif

.PRECIOUS: .build/%.o

bin/%: .build/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS)

%.so: .build/%.o
	$(CXX) $(LDFLAGS) -shared $(filter %.o,$^) -o $@ $(LDLIBS)

.build/%.o: src/%.cc
	@mkdir -pv $(dir $@)
	$(CXX) -std=$(CPPSTD) $(CFLAGS) $(DEPFLAGS) -c $(filter %.cc,$^) -o $@

.build/%.o: src/%.c
	@mkdir -pv $(dir $@)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $(filter %.c,$^) -o $@

-include $(shell [ -d '.build' ] && find .build -type f -name '*.d')

endif ###############################################################

clean:
	@rm -frv .build bin

