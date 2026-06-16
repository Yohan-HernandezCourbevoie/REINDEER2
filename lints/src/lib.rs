#![feature(rustc_private)]
#![warn(unused_extern_crates)]

// extern crate rustc_arena;
// extern crate rustc_ast;
// extern crate rustc_ast_pretty;
// extern crate rustc_data_structures;
// extern crate rustc_errors;
extern crate rustc_hir;
// extern crate rustc_hir_pretty;
// extern crate rustc_index;
// extern crate rustc_infer;
// extern crate rustc_lexer;
extern crate rustc_middle;
// extern crate rustc_mir_dataflow;
// extern crate rustc_parse;
// extern crate rustc_span;
// extern crate rustc_target;
// extern crate rustc_trait_selection;

use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;

dylint_linting::declare_late_lint! {
    /// ### What it does
    /// Warns when `seq` is called on `bio::io::fasta::Record`.
    ///
    /// ### Why is this bad?
    /// `seq` is a footgun to use because REINDEER2 will treat differently uppercase and lowercase ascii. Prefer `uppercase_seq` instead.
    ///
    /// ### Example
    ///
    /// ```rust
    /// let r: Record = /* ... */ ;
    /// let s = r.seq(); // warned
    /// ```
    ///
    /// Use instead:
    ///
    /// ```rust
    /// let r: Record = /* ... */ ;
    /// let s = r.uppercase_seq(); // fine
    /// ```
    pub RECORD_FOOTGUN,
    Warn,
    "avoid calling `seq` on `Record`"
}

impl<'tcx> LateLintPass<'tcx> for RecordFootgun {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'_>) {
        if let ExprKind::MethodCall(path, receiver, _, span) = &expr.kind {
            if path.ident.name.as_str() != "seq" {
                return;
            }

            let ty = cx.typeck_results().expr_ty(receiver).peel_refs();

            if let ty::Adt(adt_def, _) = ty.kind() {
                let def_path = cx.tcx.def_path_str(adt_def.did());
                if def_path == "bio::io::fasta::Record" {
                    cx.tcx.dcx().span_warn(
                        *span,
                        "avoid `seq` on `Record`, use `uppercase_seq()` instead",
                    );
                }
            }
        }
    }
}

#[test]
fn ui() {
    dylint_testing::ui_test(env!("CARGO_PKG_NAME"), "ui");
}
