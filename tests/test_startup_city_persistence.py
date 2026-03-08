from __future__ import annotations

import ast
from pathlib import Path


def test_save_last_city_called_only_in_startup_resolve_city() -> None:
    """Guard against reintroducing duplicate startup city persistence."""
    project_root = Path(__file__).resolve().parents[1]
    src = (project_root / "src" / "zstarview" / "startup.py").read_text(encoding="utf-8")
    tree = ast.parse(src)

    save_calls: list[tuple[str | None, int]] = []

    class Visitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.func_stack: list[str] = []

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            self.func_stack.append(node.name)
            self.generic_visit(node)
            self.func_stack.pop()

        def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
            self.func_stack.append(node.name)
            self.generic_visit(node)
            self.func_stack.pop()

        def visit_Call(self, node: ast.Call) -> None:
            if isinstance(node.func, ast.Name) and node.func.id == "save_last_city":
                current_func = self.func_stack[-1] if self.func_stack else None
                save_calls.append((current_func, node.lineno))
            self.generic_visit(node)

    Visitor().visit(tree)

    assert len(save_calls) == 1, f"Expected 1 save_last_city call, got {save_calls}"
    assert save_calls[0][0] == "_startup_resolve_city", save_calls
