use std::fs::File;
use std::io::{self, BufWriter};
use std::ops;
use io::Write;

type Color = u32;

const WIDTH: usize = 500;
const HEIGHT: usize = 250;
const BACKGROUND: Color = 0x00FF00;

#[derive(PartialEq, Clone, Copy)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64
}

impl Vec3 {
    fn new(x: f64, y: f64, z:f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    fn fromf(x: f64) -> Vec3{
        Vec3::new(x, x, x)
    }

    fn fromi(x: i64) -> Vec3{
        Vec3::fromf(x as f64)
    }

    fn zero() -> Vec3 {
        Vec3::fromi(0)
    }

    fn one() -> Vec3 {
        Vec3::fromi(1)
    }

    fn color(&self) -> Color {
        let mut res: Color = 0; 
        res += ((self.x.min(1.0).max(0.0) * 255.01) as u32) << 8 * 2;
        res += ((self.y.min(1.0).max(0.0) * 255.01) as u32) << 8 * 1;
        res += ((self.z.min(1.0).max(0.0) * 255.01) as u32) << 8 * 0;
        res
    }

    fn sum(&self) -> f64 {
        self.x + self.y + self.z
    }

    fn sqr_len(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn len(&self) -> f64 {
        self.sqr_len().sqrt()
    }

    fn norm(&self) -> Vec3 {
        *self / Vec3::fromf(self.len())
    }

    fn dot(&self, other: &Vec3) -> f64 {
        (*self * *other).sum()
    }

    fn cross(&self, other: &Vec3) -> Vec3 {
        Vec3::new(self.y * other.z - self.z * other.y,
                  -(self.x * other.z - self.z * other.x),
                  self.x * other.y - self.y * other.x)
    }


}


impl ops::Add<Vec3> for Vec3 {
    type Output = Self;

    fn add(self, other: Vec3) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

}

impl ops::AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        *self = *self * other;
    }
}

impl ops::Add<Vec3> for f64 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Self::Output {
        other + self
    }
}

impl ops::Add<f64> for Vec3 {
    type Output = Self;

    fn add(self, other: f64) -> Self::Output {
        self + Self::fromf(other)
    }
}

impl ops::AddAssign<f64> for Vec3 {
    fn add_assign(&mut self, other: f64) {
        *self = *self + other;
    }
}

impl ops::Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl ops::Sub<Vec3> for Vec3 {
    type Output = Self;
    
    fn sub(self, other: Vec3) -> Self::Output {
        self + (-other)
    }
}

impl ops::SubAssign<Vec3> for Vec3 {
    fn sub_assign(&mut self, other: Vec3) {
        *self = *self - other;
    }
}

impl ops::Sub<f64> for Vec3 {
    type Output = Self;

    fn sub(self, other: f64) -> Self::Output {
        self + (-other)
    }
}

impl ops::SubAssign<f64> for Vec3 {
    fn sub_assign(&mut self, other: f64) {
        *self = *self - other;
    }
}

impl ops::Mul<Vec3> for Vec3 {
    type Output = Self;
    
    fn mul(self, other: Vec3) -> Self::Output {
        Self::new(self.x * other.x, self.y * other.y, self.z * other.z)
    }
}

impl ops::MulAssign<Vec3> for Vec3 {
    fn mul_assign(&mut self, other: Vec3) {
        *self = *self * other;
    }
}

impl ops::Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Self::Output {
        other * self
    }

}

impl ops::Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        self * Self::fromf(other)
    }
    
}

impl ops::MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, other: f64) {
        *self = *self * other;
    }
}

impl ops::Div<Vec3> for Vec3 {
    type Output = Self;

    fn div(self, other: Vec3) -> Self::Output {
        Self::new(self.x / other.x, self.y / other.y, self.z / other.z)
    }
}

impl ops::DivAssign<Vec3> for Vec3 {
    fn div_assign(&mut self, other: Vec3) {
        *self = *self / other;
    }
}

impl ops::Div<f64> for Vec3 {
    type Output = Self;
    
    fn div(self, other: f64) -> Self::Output {
        self / Self::fromf(other)
    }
}

impl ops::DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, other: f64) {
        *self = *self / other;
    }
}

struct Ray {
    a: Vec3,
    b: Vec3
}

impl Ray {
    fn new(a: Vec3, b: Vec3) -> Ray {
        Ray {a, b}
    }

    fn origin(&self) -> Vec3 {
        self.a
    }

    fn direction(&self) -> Vec3 {
        self.b
    }

    fn point_at_t(&self, t: f64) -> Vec3 {
        self.a + self.b * t
    }
}

fn background_render(r: &Ray) -> Color {
    let norm_direction = r.direction().norm();
    let t = 0.5 * (norm_direction.y + 1.0);
    ((1.0 - t) * Vec3::one() + t * Vec3::new(0.5, 0.7, 1.0)).color()
}

fn render<T: Hitable>(ray: &Ray, hitable: &T) -> Color {
    let mut hit_rec = HitRecord::Empty();
    if hitable.hit(ray, 0.0, f64::MAX, &mut hit_rec) {
        (0.5 * (hit_rec.normal + 1.0)).color()
    }
    else {
        background_render(ray)
    }
}

#[derive(Clone, Copy)]
struct HitRecord {
    t: f64,
    p: Vec3,
    normal: Vec3
}

impl HitRecord {
    fn Empty() -> HitRecord {
        HitRecord { t: 0.0, p: Vec3::zero(), normal: Vec3::zero()}
    }
}

trait Hitable {            
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64, hit_rec: &mut HitRecord) -> bool;
}

struct Sphere {
    center: Vec3,
    radius: f64
}

impl Sphere {
    fn new(center: Vec3, radius: f64) -> Sphere{
        Sphere { center, radius }
    }
}

impl Hitable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64, hit_rec: &mut HitRecord) -> bool {
        let oc = ray.origin() - self.center;
        let a = ray.direction().dot(&ray.direction());
        let b = oc.dot(&ray.direction());
        let c = oc.dot(&oc) - self.radius * self.radius;
        let discriminant = b*b - a*c;
        if discriminant > 0.0 {
            let temp = (-b - (b*b-a*c).sqrt()) / a;
            if temp < t_max && temp > t_min {
                hit_rec.t = temp;
                hit_rec.p = ray.point_at_t(hit_rec.t);
                hit_rec.normal = (hit_rec.p - self.center) / self.radius;
                return true
            }
            let temp = (-b + (b*b-a*c).sqrt()) / a;
            if temp < t_max && temp > t_min {
                hit_rec.t = temp;
                hit_rec.p = ray.point_at_t(hit_rec.t);
                hit_rec.normal = (hit_rec.p - self.center) / self.radius;
                return true
            }
        }
        return false;
    }
    
}

struct HitableList {
    list: Vec<HitableThing>
}

impl Hitable for HitableList {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64, hit_rec: &mut HitRecord) -> bool {
        let mut tmp_rec = HitRecord::Empty();
        let mut hit_anything = false;
        let mut closest: f64 = t_max;
        for hitable_thing in self.list.iter() {
            if hitable_thing.hit(ray, t_min, closest, &mut tmp_rec) {
                hit_anything = true;
                closest = tmp_rec.t;
                *hit_rec = tmp_rec;
            }

        }
        return hit_anything;
    }
}

enum HitableThing {
    Sphere(Sphere),
    HitableList(HitableList)
}

impl Hitable for HitableThing {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64, hit_rec: &mut HitRecord) -> bool {
        match self {
            HitableThing::Sphere(sphere) => sphere.hit(ray, t_min, t_max, hit_rec),
            HitableThing::HitableList(sphere) => sphere.hit(ray, t_min, t_max, hit_rec),
        }
    }
}



fn main() {
    let mut pixels = [BACKGROUND; WIDTH * HEIGHT];
    let ll_corner = Vec3::new(-2.0, -1.0, -1.0);
    let hor = Vec3::new(4.0, 0.0, 0.0);
    let vert = Vec3::new(0.0, 2.0, 0.0);
    let origin = Vec3::zero();
    let mut hitable_list = HitableList {list: Vec::new()};
    hitable_list.list.push(
        HitableThing::Sphere(Sphere::new(Vec3::new(0.0,0.0,-1.0), 0.5))
    );
    hitable_list.list.push(
        HitableThing::Sphere(Sphere::new(Vec3::new(0.0,-100.5,-1.0), 100.0))
    );
    for i in (0..HEIGHT).rev() {
        for j in 0..WIDTH {
            let u = j as f64 / WIDTH as f64;
            let v = i as f64 / HEIGHT as f64;
            let ray = Ray::new(origin, ll_corner + u*hor + v*vert);
            let col = render(&ray, &hitable_list);
            // println!("{:#08x}", col);
            pixels[(HEIGHT - i - 1) * WIDTH + j] = col;
        }
    }
    save_as_ppm("test.ppm", &pixels);
}

fn write_all<W: Write, I: Iterator<Item=u8>>(writer: &mut W, iter: I) -> io::Result<()> {
    let mut buffer = BufWriter::with_capacity(1024, writer);

    for el in iter {
        buffer.write(&[el])?;
    }
    Ok(())
}

fn save_as_ppm(file_path: &str, pixels: &[Color]) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    write!(file, "P6\n{} {} 255\n", WIDTH, HEIGHT)?;
    let to_write = (0..WIDTH * HEIGHT * 3).map(|i| {
        let pixel = pixels[i / 3];
        ((pixel >> 8 * (2 - i % 3))) as u8
    });
    write_all(&mut file, to_write)?;
    println!("Saved {}", file_path);
    Ok(())
}