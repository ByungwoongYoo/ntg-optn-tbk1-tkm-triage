import Lake

open Lake DSL

package U4Formal where
  version := v!"2.0.0"

@[default_target]
lean_lib U4Formal where
  roots := #[`U4Formal]
