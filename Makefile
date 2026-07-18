ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

BUILD_DIR ?= $(ROOT_DIR)/build
INSTALL_PREFIX ?=
MROCK_EXTRA_INCLUDE_DIRS ?=
CMAKE_PREFIX_PATH ?=

CONFIGURE_STAMP := $(BUILD_DIR)/CMakeCache.txt

define CMAKE_CONFIGURE_CMD
cmake -S $(ROOT_DIR) -B $(BUILD_DIR) \
	$(if $(INSTALL_PREFIX),-DCMAKE_INSTALL_PREFIX="$(INSTALL_PREFIX)") \
	$(if $(MROCK_EXTRA_INCLUDE_DIRS),-DMROCK_EXTRA_INCLUDE_DIRS="$(MROCK_EXTRA_INCLUDE_DIRS)") \
	$(if $(CMAKE_PREFIX_PATH),-DCMAKE_PREFIX_PATH="$(CMAKE_PREFIX_PATH)")
endef

$(CONFIGURE_STAMP):
	mkdir -p $(BUILD_DIR)
	$(CMAKE_CONFIGURE_CMD)

configure:
	mkdir -p $(BUILD_DIR)
	$(CMAKE_CONFIGURE_CMD)

reconfigure:
	rm -f $(CONFIGURE_STAMP)
	$(MAKE) configure \
		BUILD_DIR="$(BUILD_DIR)" \
		INSTALL_PREFIX="$(INSTALL_PREFIX)" \
		MROCK_EXTRA_INCLUDE_DIRS="$(MROCK_EXTRA_INCLUDE_DIRS)" \
		CMAKE_PREFIX_PATH="$(CMAKE_PREFIX_PATH)"

build: $(CONFIGURE_STAMP)
	cmake --build $(BUILD_DIR) --parallel

test: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure

install: build
	cmake --install $(BUILD_DIR) \
		$(if $(INSTALL_PREFIX),--prefix "$(INSTALL_PREFIX)")

all: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure
	cmake --install $(BUILD_DIR) \
		$(if $(INSTALL_PREFIX),--prefix "$(INSTALL_PREFIX)")

clean:
	rm -rf $(BUILD_DIR) $(INSTALL_PREFIX)

.PHONY: all build test install clean configure reconfigure