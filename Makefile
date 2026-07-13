ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD_DIR := $(ROOT_DIR)/build
INSTALL_PREFIX := $(ROOT_DIR)/install
CONFIGURE_STAMP := $(BUILD_DIR)/CMakeCache.txt

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(CONFIGURE_STAMP):
	mkdir -p $(BUILD_DIR)
	cmake -S $(ROOT_DIR) -B $(BUILD_DIR)

test: $(CONFIGURE_STAMP)
	cmake --build $(BUILD_DIR) --parallel
	ctest --test-dir $(BUILD_DIR) --output-on-failure

install: $(CONFIGURE_STAMP)
	@cmake --install $(BUILD_DIR)

all:
	@bash $(ROOT_DIR)/scripts/build_and_test.sh

clean:
	@rm -rf $(BUILD_DIR) $(INSTALL_PREFIX)

.PHONY: all test install clean
