// Thankfully provided by github.com/jorson
// https://github.com/typst/typst/issues/380#issuecomment-1523884719
// https://typst.app/project/roACTA6jaO8lqADEG5PFc-


#let foldl1(a, f) = a.slice(1).fold(a.first(), f)
#let concat(a) = foldl1(a, (acc, x) => acc + x)
#let nonumber(e) = math.equation(block: true, numbering: none, e)

#let eq(es, numberlast: false ) = if es.has("children") {
  let esf = es.children.filter(x => not ([ ], [#parbreak()]).contains(x))
  let bodyOrChildren(e) = if e.body.has("children") { concat(e.body.children) } else { e.body }
  let hideEquation(e) = if e.has("numbering") and e.numbering == none {
    nonumber(hide(e))
  } else [
    $ #hide(bodyOrChildren(e)) $ #{if e.has("label") { e.label }}
  ]
  let hidden = box(concat(
    if numberlast == true {
      esf.slice(0, esf.len()-1).map(e => nonumber(hide(e))) + (hideEquation(esf.last()),)
    } else if numberlast == false {
      esf.map(e => hideEquation(e))
    } else if numberlast == none {
      esf.map(e => nonumber(hide(e)))
    }))
  let folder(acc, e) = acc + if acc != [] { linebreak() } + e
  let aligned = math.equation(block: true, numbering: none, esf.fold([], folder))

  hidden
  // style(s => v(-measure(hidden, s).height, weak: true))
  context {
    v(-measure(hidden).height, weak: true)
  }
  aligned
}
