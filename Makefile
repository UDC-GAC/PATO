BUILD=release

INCLUDES+=-I include

CXX=g++
CXXFLAGS+=-std=c++14

cxxflags.release=-O3 -march=native -DNDEBUG

CXXFLAGS+=$(cxxflags.$(BUILD))

LD=g++

ldflags.release=

LDFLAGS+=$(ldflags.$(BUILD))

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
