NAME = CurvatureEngine

CC = clang++

CFLAGS = -g -std=c++17 -O3\
         -mavx2 -mfma -march=native -mtune=native \
         -funroll-loops -fvectorize -ffp-contract=fast \
         -freciprocal-math -ffast-math -fstrict-aliasing \
         -fomit-frame-pointer -flto=full -mprefer-vector-width=256 -fopenmp \
         -I$(INC_DIR)

SRC_DIR = srcs
INC_DIR = includes
OBJ_DIR = build
OUT_VTK_DIR = output

SRC = $(shell find $(SRC_DIR) -name '*.cpp')

OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))

TOTAL := $(words $(SRC))

GREEN := \033[0;32m
YELLOW := \033[0;33m
RED := \033[0;31m
NC := \033[0m

.PHONY: all clean fclean re


all: $(NAME)
	@if [ -f $(OBJ_DIR)/.counter ]; then rm $(OBJ_DIR)/.counter; fi
	@echo -e "$(GREEN)Build complete!$(NC)"

$(NAME): $(OBJ)
	@echo -e "$(YELLOW)Linking $@...$(NC)"
	$(CC) $(CFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(dir $@) 
	@if [ ! -f $(OBJ_DIR)/.counter ]; then echo 0 > $(OBJ_DIR)/.counter; fi; \
	cnt=$$(cat $(OBJ_DIR)/.counter); \
	cnt=$$((cnt+1)); \
	echo $$cnt > $(OBJ_DIR)/.counter; \
	perc=$$((100 * cnt / $(TOTAL))); \
	printf "$(GREEN)%3d%% Compiling $< [$$cnt/$(TOTAL)]$(NC)\n" $$perc; \
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	@echo -e "$(RED)Cleaning up...$(NC)"
	rm -rf $(OBJ_DIR)

fclean: clean
	rm -f $(NAME)

re: fclean all

