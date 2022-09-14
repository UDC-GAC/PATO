BUILD=release

INCLUDES+=-Iinclude

CXX=g++
CXXFLAGS+=-std=c++17 -Wall

cxxflags.serial=-flto -O3 -march=native -s
cxxflags.profile=-flto -O3 -march=native -p -g
cxxflags.release=-fopenmp -flto -O3 -march=native -s

CXXFLAGS+=$(cxxflags.$(BUILD)) -DNDEBUG -DSEQAN_DISABLE_VERSION_CHECK -DSEQAN_ENABLE_PARALLELISM=0

LD=g++

ldflags.serial=-flto
ldflags.profile=-flto -p -g
ldflags.release=-fopenmp -flto

# see file seqan/allocator_interface.h:211
LDFLAGS+=-Wno-alloc-size-larger-than $(ldflags.$(BUILD))

SRCSDIR=src
OBJSDIR=obj
DESTDIR=target

SRCS=$(wildcard $(SRCSDIR)/*.cpp)
OBJS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/$(BUILD)/%.o,$(SRCS))
DEPS=$(patsubst $(SRCSDIR)/%.cpp,$(OBJSDIR)/$(BUILD)/%.d,$(SRCS))

pato.$(BUILD): $(OBJS)
	@mkdir -p $(DESTDIR)
	$(LD) $(LDFLAGS) -o $(DESTDIR)/$@ $^ $(LDLIBS)

-include $(DEPS)

$(OBJSDIR)/$(BUILD)/%.o: $(SRCSDIR)/%.cpp
	@mkdir -p $(OBJSDIR)/$(BUILD)
	$(CXX) $(CXXFLAGS) -MMD -MP $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJSDIR) $(DESTDIR)
