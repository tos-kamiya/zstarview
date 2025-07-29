import os

import pygame

_dir = os.path.dirname(os.path.abspath(__file__))
print(f"{_dir=}")

font_path = os.path.join(_dir, "Noto_Sans_Symbols", "NotoSansSymbols-VariableFont_wght.ttf")

pygame.font.init()
font = pygame.font.Font(font_path, 24)
text = font.render("☉☿♀♂♃♄12345", True, (255, 255, 0))
pygame.image.save(text, "test_text.png")
