ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

BUILD_DIR ?= $(ROOT_DIR)/build
INSTALL_PREFIX ?=

CONFIGURE_STAMP := $(BUILD_DIR)/CMakeCache.txt

$(CONFIGURE_STAMP):
	mkdir -p $(BUILD_DIR)
	cmake -S $(ROOT_DIR) -B $(BUILD_DIR) \
	    $(if $(INSTALL_PREFIX),-DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX))

build: $(CONFIGURE_STAMP)
	cmake --build $(BUILD_DIR) --parallel

test: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure

install: build
	cmake --install $(BUILD_DIR) \
	    $(if $(INSTALL_PREFIX),--prefix $(INSTALL_PREFIX))

all: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure
	cmake --install $(BUILD_DIR) \
	    $(if $(INSTALL_PREFIX),--prefix $(INSTALL_PREFIX))

clean:
	rm -rf $(BUILD_DIR) $(INSTALL_PREFIX)

.PHONY: all build test install clean