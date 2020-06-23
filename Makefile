TARGET = romberg-par

SRCDIR = src
SOURCE = $(wildcard $(SRCDIR)/*.cpp)
INCDIR = $(SRCDIR)/include

OBJECTS := $(patsubst %.cpp,%.o, $(SOURCE))

CPPFLAGS = -O3
CXXFLAGS = $(CPPFLAGS) -std=c++11 -I$(INCDIR)
LDFLAGS = -lpthread

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	@$(CXX) $^ -o $@ $(LDFLAGS)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp*
	@echo "Compiling $<.."
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	@rm $(TARGET)

.PHONY: clean