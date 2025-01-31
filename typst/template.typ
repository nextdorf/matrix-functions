// Based on preview/arkheion:0.1.0, https://github.com/mgoulao/arkheion/tree/main, https://typst.app/universe/package/arkheion/, MIT liensed

#let minimal-document(
  title: "",
  abstract: [],
  keywords: (),
  authors: (),
  date: none,
  body,
) = {
  // Set the document's basic properties.
  set document(author: authors.map(a => a.name), title: title)
  set page(
    margin: (left: 12.5mm, right: 12.5mm, top: 12.5mm, bottom: 12.5mm),
    numbering: "1",
    number-align: center,
  )
  set text(font: "New Computer Modern", lang: "en")
  show math.equation: set text(weight: 400)
  show math.equation: set block(spacing: 0.65em)
  set math.equation(
    numbering: "(1)",
    supplement: none,
    )
  set heading(numbering: "1.1")

  let refcol = color.hsl(245deg, 80%, 45%)
  show link: set text(refcol)
  show ref: it => {
    set text(refcol)
    if it.element != none and it.element.func() == math.equation {
      link(it.target)[(#it)]
    } else {
      it
    }
  }
  // Set run-in subheadings, starting at level 4.
  show heading: it => {
    // H1 and H2
    if it.level == 1 {
      pad(
        top: .5em,
        // bottom: 1em,
        it
      )
    }
    else if it.level == 2 {
      pad(
        top: .5em,
        bottom: .25em,
        it
      )
    }
    else if it.level > 3 {
      text(11pt, weight: "bold", it.body + " ")
    } else {
      it
    }
  }

  line(length: 100%, stroke: 2pt)
  // Title row.
  pad(
    bottom: 4pt,
    top: 4pt,
    align(center)[
      #block(text(weight: 500, 1.75em, title))
      #v(1em, weak: true)
    ]
  )
  line(length: 100%, stroke: 2pt)

  // Author information.
  pad(
    top: 0.5em,
    x: 2em,
    grid(
      columns: (1fr,) * calc.min(3, authors.len()),
      gutter: 1em,
      ..authors.map(author => align(center)[
        *#author.name* \
        #if author.keys().contains("subtitle") {
          author.subtitle
        }
      ]),
    ),
  )

  // align(center)[#date]

  // Table of Contents.
  outline()

  // Abstract.
  if abstract != [] {    
    pad(
      x: 3em,
      top: 1em,
      bottom: 0.4em,
      align(center)[
        #heading(
          outlined: false,
          numbering: none,
          text(0.85em, smallcaps[Abstract]),
        )
        #set par(justify: true)
        #set text(hyphenate: false)

        #abstract
      ],
    )
  }

  // // Keywords
  // if keywords.len() > 0 {
  //     [*_Keywords_* #h(0.3cm)] + keywords.map(str).join(" · ")
  // }
  // Main body.
  set par(justify: true)
  set text(hyphenate: false)

  body
}

// #let arkheion-appendices(body) = {
//   counter(heading).update(0)
//   counter("appendices").update(1)

//   set heading(
//     numbering: (..nums) => {
//       let vals = nums.pos()
//       let value = "ABCDEFGHIJ".at(vals.at(0) - 1)
//       if vals.len() == 1 {
//         return "APPENDIX " + value
//       }
//       else {
//         return value + "." + nums.pos().slice(1).map(str).join(".")
//       }
//     }
//   );
//   [#pagebreak() #body]
// }