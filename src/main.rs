fn main() {
    let xmin : f64 = 0.;
    let xmax : f64 = 2.;
    let ymin : f64 = 0.;
    let ymax : f64 = 2.;

    let nx : f64 = 21.;  // These should be unsigned ints, but rust is
    let ny : f64 = 21.;  // toooo explicit

    let dx = ( xmax - xmin ) / nx;
    let dy = ( ymax - ymin ) / ny;

    println!("dx = {}\ndy = {}\n", dx, dy);

    let mut u : [f64; nx];
}
