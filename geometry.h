#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

const double delta = 1e-5;

bool compare(double x, double y) {
  return (abs(x - y) < delta);
}

struct Point {
    double x_;
    double y_;

    Point(double x, double y) : x_(x), y_(y) {}

    Point() : x_(0), y_(0) {}
};

bool operator==(const Point &point_1, const Point &point_2) {
  return (compare(point_1.x_, point_2.x_) && compare(point_1.y_, point_2.y_));
}

bool operator!=(const Point &point_1, const Point &point_2) {
  return !(point_1 == point_2);
}

double dist(const Point &point_1, const Point &point_2) {
  double ans = sqrt(pow(point_1.x_ - point_2.x_, 2) + pow(point_1.y_ - point_2.y_, 2));
  return ans;
}

Point mid(const Point &pount_1, const Point &point_2) {
  Point ans((pount_1.x_ + point_2.x_) / 2, (pount_1.y_ + point_2.y_) / 2);
  return ans;
}


class Line {
public:
    Line(const Point &point_1, const Point &point_2) {
      a_ = point_2.y_ - point_1.y_;
      b_ = point_1.x_ - point_2.x_;
      c_ = point_1.y_ * point_2.x_ - point_1.x_ * point_2.y_;
    }

    Line() : a_(0), b_(0), c_(0) {}

    Line(double k, double b_1) : a_(k), b_(-1), c_(b_1) {
    }

    Line(const Point &p, double k) : a_(k), b_(-1), c_(p.y_ - k * p.x_) {
    }

    Line(double a, double b, double c) : a_(a), b_(b), c_(c) {}

    double get_a() const {
      return a_;
    }

    double get_b() const {
      return b_;
    }

    double get_c() const {
      return c_;
    }

private:
    double a_;
    double b_;
    double c_;
};

bool operator==(const Line &line_1, const Line &line_2) {
  double a_1 = line_1.get_a();
  double a_2 = line_2.get_a();
  double b_1 = line_1.get_b();
  double b_2 = line_2.get_b();
  double c_1 = line_1.get_c();
  double c_2 = line_2.get_c();
  if (a_1 != 0) {
    double d = a_2 / a_1;
    return (compare(b_2, d * b_1) && compare(c_2, d * c_1);
  }
  double d = b_2 / b_1;
  return (compare(a_2, d * a_1) && compare(c_2, d * c_1));
}

bool operator!=(const Line &line_1, const Line &line_2) {
  return !(line_1 == line_2);
}

Point crossing(const Line &line_1, const Line &line_2) {
  Point ans;
  if (line_2.get_a() == 0) {
    ans.y_ = -line_2.get_c() / line_2.get_b();
    ans.x_ = (-line_1.get_c() - line_1.get_b() * ans.y_) / line_1.get_a();
    return ans;
  }
  if (line_1.get_a() == 0) {
    ans.y_ = -line_1.get_c() / line_1.get_b();
    ans.x_ = (-line_2.get_c() - line_2.get_b() * ans.y_) / line_2.get_a();
    return ans;
  }
  double d = line_2.get_a() / line_1.get_a();
  ans.y_ = (-line_1.get_c() * d + line_2.get_c()) / (d * line_1.get_b() - line_2.get_b());
  ans.x_ = (-line_2.get_c() - line_2.get_b() * ans.y_) / line_2.get_a();
  return ans;
}

Line normal(const Line &line_1, const Point &point) {
  if (line_1.get_b() == 0) {
    Line L_2(0, 1, -point.y_);
    return L_2;
  }
  if (line_1.get_a() == 0) {
    Line L_2(1, 0, -point.x_);
    return L_2;
  }
  Line L_2(point, line_1.get_b() / line_1.get_a());
  return L_2;
}

Line bisectrix(const Line &line_1, const Line &line_2) {
  double z_1 = sqrt(pow(line_1.get_a(), 2) + pow(line_1.get_b(), 2));
  double z_2 = sqrt(pow(line_2.get_a(), 2) + pow(line_2.get_b(), 2));
  Line ans(line_1.get_a() * z_2 - line_2.get_a() * z_1, line_1.get_b() * z_2 - line_2.get_b() * z_1,
           line_1.get_c() * z_2 - line_2.get_c() * z_1);
  return ans;
}


class Shape {
public:
    virtual double perimeter() = 0;

    virtual double area() = 0;

    virtual void rotate(Point center, double angle) = 0;

    virtual void reflex(Point center) = 0;

    virtual void reflex(Line axis) = 0;

    virtual void scale(Point center, double coefficient) = 0;

    virtual bool containsPoint(Point point) = 0;

    virtual bool isSimilarTo(const Shape &another) = 0;

    virtual bool operator==(const Shape &another) const = 0;

    virtual bool operator!=(const Shape &another) const = 0;

    virtual bool isCongruentTo(const Shape &another) = 0;

    virtual ~Shape() = 0;
};

Shape::~Shape() {}

class Polygon : public Shape {
public:
    int verticesCount() const {
      return vertices_.size();
    }

    std::vector<Point> getVertices() const {
      std::vector<Point> ans;
      for (int i = 0; i < int(vertices_.size()); ++i) {
        ans.push_back(vertices_[i]);
      }
      return ans;
    };

    bool isConvex() {
      Point a;
      a.x_ = vertices_[vertices_.size() - 1].x_ - vertices_[0].x_;
      a.y_ = vertices_[vertices_.size() - 1].y_ - vertices_[0].y_;
      Point b;
      b.x_ = vertices_[0].x_ - vertices_[1].x_;
      b.y_ = vertices_[0].y_ - vertices_[1].y_;
      double canon = a.x_ * b.y_ - a.y_ * b.x_;
      for (int i = 1; i < int(vertices_.size()); ++i) {
        a.x_ = vertices_[i - 1].x_ - vertices_[i].x_;
        a.y_ = vertices_[i - 1].y_ - vertices_[i].y_;
        b.x_ = vertices_[i].x_ - vertices_[(i + 1) % vertices_.size()].x_;
        b.y_ = vertices_[i].y_ - vertices_[(i + 1) % vertices_.size()].y_;
        if ((a.x_ * b.y_ - a.y_ * b.x_) * canon < 0) {
          return false;
        }
      }
      return true;
    }

    double perimeter() override {
      double ans = dist(vertices_[vertices_.size() - 1], vertices_[0]);
      for (int i = 1; i < int(vertices_.size()); ++i) {
        ans += dist(vertices_[i - 1], vertices_[i]);
      }
      return ans;
    }

    double area() override {
      double ans = 0;
      for (int i = 0; i < int(vertices_.size()); ++i) {
        Point point_1 = vertices_[(i + 1) % vertices_.size()];
        Point point_2 = vertices_[i];
        ans += (point_1.x_ - point_2.x_) * (point_1.y_ + point_2.y_);
      }
      return ans / 2;
    }

    void rotate(Point center, double angle) override {
      double c = cos(angle);
      double s = sin(angle);
      for (int i = 0; i < int(vertices_.size()); ++i) {
        Point F_0 = vertices_[i];
        vertices_[i].x_ = center.x_ + (F_0.x_ - center.x_) * c - (F_0.y_ - center.y_) * s;
        vertices_[i].y_ = center.x_ + (F_0.x_ - center.x_) * s + (F_0.y_ - center.y_) * c;
      }
    }

    void reflex(Point center) override {
      for (int i = 0; i < int(vertices_.size()); ++i) {
        Point point = vertices_[i];
        vertices_[i].x_ = 2 * center.x_ - point.x_;
        vertices_[i].y_ = 2 * center.y_ - point.y_;
      }
    }

    void reflex(Line axis) override {
      for (int i = 0; i < int(vertices_.size()); ++i) {
        Line line_1 = normal(axis, vertices_[i]);
        Point center = crossing(line_1, axis);
        Point point = vertices_[i];
        vertices_[i].x_ = 2 * center.x_ - point.x_;
        vertices_[i].y_ = 2 * center.y_ - point.y_;
      }
    }

    void scale(Point center, double coefficient) override {
      for (int i = 0; i < int(vertices_.size()); ++i) {
        Point point = vertices_[i];
        vertices_[i].x_ = center.x_ + coefficient * (point.x_ - center.x_);
        vertices_[i].y_ = center.y_ + coefficient * (point.y_ - center.y_);
      }
    }

    virtual bool containsPoint(Point point) override {
      int size = vertices_.size() - 1;
      int count = 0;
      for (int i = 0; i < int(vertices_.size()); ++i) {
        double x_1 = vertices_[i].x_;
        double x_2 = vertices_[size].x_;
        double x_3 = point.x_;
        double y_1 = vertices_[i].y_;
        double y_2 = vertices_[size].y_;
        double y_3 = point.y_;
        if ((y_1 < y_3 && y_2 >= y_3) || (y_1 >= y_3 && y_2 < y_3)) {
          if (x_1 + (y_3 - y_1) / (y_2 - y_1) * (x_2 - x_1) < x_3) {
            ++count;
          }
        }
        size = i;
      }
      return count % 2;
    }

    bool operator==(const Shape &another) const override {
      const Polygon *p = dynamic_cast<const Polygon *>(&another);
      if (p) {
        if (vertices_.size() != p->vertices_.size()) {
          return false;
        }
        int s = vertices_.size();
        for (int i = 0; i < s; ++i) {
          int j = 0;
          while (j < s) {
            if (vertices_[(i + j) % s] == p->vertices_[j]) {
              ++j;
            } else {
              break;
            }
            if (j == s) {
              return true;
            }
          }
        }
        for (int i = 0; i < s; ++i) {
          int j = s - 1;
          int t = 0;
          while (t < s) {
            if (vertices_[(i + t) % s] == p->vertices_[j]) {
              --j;
              ++t;
            } else {
              break;
            }
            if (t == s) {
              return true;
            }
          }
        }
      }
      return false;
    }

    bool operator!=(const Shape &another) const override {
      return !(*this == another);
    }


    bool isCongruentTo(const Shape &another) override {
      const Polygon *p = dynamic_cast<const Polygon *>(&another);
      if (p) {
        if (vertices_.size() != p->vertices_.size()) {
          return false;
        }
        int s = vertices_.size();
        for (int i = 0; i < s; ++i) {
          int j = 0;
          while (j < s) {
            Point point_1 = vertices_[(i + j) % s];
            Point point_2 = vertices_[(i + j + 1) % s];
            Point point_3 = vertices_[(i + j + 2) % s];
            Point point_4 = p->vertices_[j];
            Point point_5 = p->vertices_[(j + 1) % s];
            Point point_6 = p->vertices_[(j + 2) % s];
            if (compare(dist(point_1, point_2), dist(point_4, point_5)) &&
                compare(dist(point_2, point_3), dist(point_5, point_6))
                && compare(dist(point_1, point_3), dist(point_4, point_6))) {
              ++j;
            } else {
              break;
            }
            if (j == s) {
              return true;
            }
          }
        }
        for (int i = 0; i < s; ++i) {
          int t = 0;
          int j = s - 1;
          while (t < s) {
            Point point_1 = vertices_[(i + t) % s];
            Point point_2 = vertices_[(i + t + 1) % s];
            Point point_3 = vertices_[(i + t + 2) % s];
            Point point_4 = p->vertices_[j % s];
            Point point_5 = p->vertices_[(j - 1 + s) % s];
            Point point_6 = p->vertices_[(j - 2 + s) % s];
            if (compare(dist(point_1, point_2), dist(point_4, point_5)) &&
                compare(dist(point_2, point_3), dist(point_5, point_6))
                && compare(dist(point_1, point_3), dist(point_4, point_6))) {
              --j;
              ++t;
            } else {
              break;
            }
            if (t == s) {
              return true;
            }
          }
        }
      }
      return false;
    }

    bool isSimilarTo(const Shape &another) override {
      const Polygon *p = dynamic_cast<const Polygon *>(&another);
      if (p) {
        if (vertices_.size() != p->vertices_.size()) {
          return false;
        }
        int s = vertices_.size();
        for (int i = 0; i < s; ++i) {
          double canon = dist(vertices_[i], vertices_[(i + 1) % s]) / dist(p->vertices_[0], p->vertices_[1]);
          int j = 0;
          while (j < s) {
            Point point_1 = vertices_[(i + j) % s];
            Point point_2 = vertices_[(i + j + 1) % s];
            Point point_3 = vertices_[(i + j + 2) % s];
            Point point_4 = p->vertices_[j];
            Point point_5 = p->vertices_[(j + 1) % s];
            Point point_6 = p->vertices_[(j + 2) % s];
            if (compare(dist(point_1, point_2) / dist(point_4, point_5), canon) &&
                compare(dist(point_2, point_3) / dist(point_5, point_6), canon) &&
                compare(dist(point_1, point_3) / dist(point_4, point_6), canon)) {
              ++j;
            } else {
              break;
            }
            if (j == s) {
              return true;
            }
          }
        }
        for (int i = 0; i < s; ++i) {
          double canon = dist(vertices_[i], vertices_[(i + 1) % s]) / dist(p->vertices_[s - 1], p->vertices_[s - 2]);
          int t = 0;
          int j = s - 1;
          while (t < s) {
            Point point_1 = vertices_[(i + t) % s];
            Point point_2 = vertices_[(i + t + 1) % s];
            Point point_3 = vertices_[(i + t + 2) % s];
            Point point_4 = p->vertices_[j % s];
            Point point_5 = p->vertices_[(j - 1 + s) % s];
            Point point_6 = p->vertices_[(j - 2 + s) % s];
            if (compare(dist(point_1, point_2) / dist(point_4, point_5), canon) &&
                compare(dist(point_2, point_3) / dist(point_5, point_6), canon) &&
                compare(dist(point_1, point_3) / dist(point_4, point_6), canon)) {
              --j;
              ++t;
            } else {
              break;
            }
            if (t == s) {
              return true;
            }
          }
        }
      }
      return false;
    }


    Polygon() {}

    Polygon(std::vector<Point> P_1) {
      for (int i = 0; i < int(P_1.size()); ++i) {
        vertices_.push_back(P_1[i]);
      }
    }

protected:
    std::vector<Point> vertices_;
};

class Ellipse : public Shape {
public:
    std::pair<Point, Point> focuses() {
      std::pair<Point, Point> ans;
      ans.first = focus_1_;
      ans.second = focus_2_;
      return ans;
    }

    double small() {
      double dis = dist(focus_1_, focus_2_) / 2;
      return sqrt(pow(d_, 2) - pow(dis, 2));
    }

    std::pair<Line, Line> directrices() {
      std::pair<Line, Line> ans;
      Line line_1(focus_1_, focus_2_);
      Line line_2 = normal(line_1, mid(focus_1_, focus_2_));
      double a = line_2.get_a();
      double b = line_2.get_b();
      double delta = d_ * (sqrt(pow(a, 2) + pow(b, 2))) / eccentricity();
      ans.first = Line(a, b, -delta);
      ans.second = Line(a, b, delta);
      return ans;
    }

    double eccentricity() {
      double dis = dist(focus_1_, focus_2_) / 2;
      double e = dis / d_;
      return e;
    }

    Point center() {
      return mid(focus_1_, focus_2_);
    }

    void rotate(Point center, double angle) override {
      double c = cos(angle);
      double s = sin(angle);
      Point focus = focus_1_;
      focus_1_.x_ = center.x_ + (focus.x_ - center.x_) * c - (focus.y_ - center.y_) * s;
      focus_1_.y_ = center.x_ + (focus.x_ - center.x_) * s + (focus.y_ - center.y_) * c;
      focus = focus_2_;
      focus_2_.x_ = center.x_ + (focus.x_ - center.x_) * c - (focus.y_ - center.y_) * s;
      focus_2_.y_ = center.x_ + (focus.x_ - center.x_) * s + (focus.y_ - center.y_) * c;
    }

    void reflex(Line axis) override {
      Line line_1 = normal(axis, focus_1_);
      Point center = crossing(line_1, axis);
      Point focus = focus_1_;
      focus_1_.x_ = 2 * center.x_ - focus.x_;
      focus_1_.y_ = 2 * center.y_ - focus.y_;
      Line line_2 = normal(axis, focus_2_);
      center = crossing(line_2, axis);
      focus = focus_2_;
      focus_2_.x_ = 2 * center.x_ - focus.x_;
      focus_2_.y_ = 2 * center.y_ - focus.y_;
    }

    void reflex(Point center) override {
      Point focus = focus_1_;
      focus_1_.x_ = 2 * center.x_ - focus.x_;
      focus_1_.y_ = 2 * center.y_ - focus.y_;
      focus = focus_2_;
      focus_2_.x_ = 2 * center.x_ - focus.x_;
      focus_2_.y_ = 2 * center.y_ - focus.y_;
    }

    void scale(Point center, double coefficient) override {
      Point focus = focus_1_;
      focus_1_.x_ = center.x_ + coefficient * (focus.x_ - center.x_);
      focus_1_.y_ = center.y_ + coefficient * (focus.y_ - center.y_);
      focus = focus_2_;
      focus_2_.x_ = center.x_ + coefficient * (focus.x_ - center.x_);
      focus_2_.y_ = center.y_ + coefficient * (focus.y_ - center.y_);
      d_ *= coefficient;
    }

    bool containsPoint(Point point) override {
      return (dist(focus_1_, point) + dist(focus_2_, point) <= 2 * d_);
    }

    double perimeter() override {
      double p = M_PI * (3 * (d_ + small()) - sqrt((3 * d_ + small()) * (d_ + 3 * small())));
      return p;
    }

    double area() override {
      return M_PI * small() * d_;
    }

    bool operator==(const Shape &another) const override {
      const Ellipse *p = dynamic_cast<const Ellipse *>(&another);
      if (p) {
        if ((focus_1_ == p->focus_1_ && focus_2_ == p->focus_2_) ||
            (focus_1_ == p->focus_2_ && focus_2_ == p->focus_1_)) {
          if (d_ == p->d_) {
            return true;
          }
        }
      }
      return false;
    }

    bool operator!=(const Shape &another) const override {
      return !(*this == another);
    }

    bool isSimilarTo(const Shape &another) override {
      const Ellipse *p = dynamic_cast<const Ellipse *>(&another);
      if (p) {
        if (d_ == p->d_ && dist(focus_1_, focus_2_) == dist(p->focus_1_, p->focus_2_)) {
          return true;
        }
      }
      return false;
    }

    bool isCongruentTo(const Shape &another) override {
      const Ellipse *p = dynamic_cast<const Ellipse *>(&another);
      if (p) {
        double delta = d_ / p->d_;
        if (dist(focus_1_, focus_2_) == delta * dist(p->focus_1_, p->focus_2_)) {
          return true;
        }
      }
      return false;
    }

    Ellipse(Point point_1, Point point_2, double d) : focus_1_(point_1), focus_2_(point_2), d_(d / 2) {}

    Ellipse() {}

protected:
    Point focus_1_;
    Point focus_2_;
    double d_;
};

class Circle : public Ellipse {
public:
    double radius() {
      return d_;
    }

    Circle(Point c, double r) {
      focus_1_ = c;
      focus_2_ = c;
      d_ = r;
    }
};

class Rectangle : public Polygon {
public:
    Point center() {
      return mid(vertices_[0], vertices_[2]);
    }

    std::pair<Line, Line> diagonals() {
      std::pair<Line, Line> ans;
      ans.first = Line(vertices_[0], vertices_[2]);
      ans.second = Line(vertices_[1], vertices_[3]);
      return ans;
    }

    Rectangle() {}

    Rectangle(Point point_1, Point point_2, double angle) {
      if (angle < 1) {
        angle = 1 / angle;
      }
      Point point_3;
      double c = cos(atan(angle));
      double s = sin(atan(angle));
      point_3.x_ = point_1.x_ + (point_2.x_ - point_1.x_) * c - (point_2.y_ - point_1.y_) * s;
      point_3.y_ = point_1.y_ + (point_2.x_ - point_1.x_) * s + (point_2.y_ - point_1.y_) * c;
      Line line_1(point_1, point_3);
      angle = 1 / angle;
      double c_1 = cos(-atan(angle));
      double s_1 = sin(-atan(angle));
      Point point_4;
      point_4.x_ = point_2.x_ + (point_1.x_ - point_2.x_) * c_1 - (point_1.y_ - point_2.y_) * s_1;
      point_4.y_ = point_2.y_ + (point_1.x_ - point_2.x_) * s_1 + (point_1.y_ - point_2.y_) * c_1;
      Line line_2(point_2, point_4);
      Point point_5 = crossing(line_1, line_2);
      Point point_6 = mid(point_1, point_2);
      Point point_7;
      point_7.x_ = 2 * point_6.x_ - point_5.x_;
      point_7.y_ = 2 * point_6.y_ - point_5.y_;
      vertices_.push_back(point_1);
      vertices_.push_back(point_5);
      vertices_.push_back(point_2);
      vertices_.push_back(point_7);
    }

};

class Square : public Rectangle {
public:
    Circle circumscribedCircle() {
      Point c = center();
      double d = dist(vertices_[0], vertices_[2]);
      Circle ans(c, d / 2);
      return ans;
    }

    Circle inscribedCircle() {
      Point c = center();
      double d = dist(vertices_[0], vertices_[1]);
      Circle ans(c, d / 2);
      return ans;
    }

    Square(const Point &point_1, const Point &point_2) {
      Point point_3;
      double angle = 1;
      double c = cos(atan(angle));
      double s = sin(atan(angle));
      point_3.x_ = point_1.x_ + (point_2.x_ - point_1.x_) * c - (point_2.y_ - point_1.y_) * s;
      point_3.y_ = point_1.y_ + (point_2.x_ - point_1.x_) * s + (point_2.y_ - point_1.y_) * c;
      Line line_1(point_1, point_3);
      double c_1 = cos(-atan(angle));
      double s_1 = sin(-atan(angle));
      Point point_4;
      point_4.x_ = point_2.x_ + (point_1.x_ - point_2.x_) * c_1 - (point_1.y_ - point_2.y_) * s_1;
      point_4.y_ = point_2.y_ + (point_1.x_ - point_2.x_) * s_1 + (point_1.y_ - point_2.y_) * c_1;
      Line line_2(point_2, point_4);
      Point point_5 = crossing(line_1, line_2);
      Point point_6 = mid(point_1, point_2);
      Point point_7;
      point_7.x_ = 2 * point_6.x_ - point_5.x_;
      point_7.y_ = 2 * point_6.y_ - point_5.y_;
      vertices_.push_back(point_1);
      vertices_.push_back(point_5);
      vertices_.push_back(point_2);
      vertices_.push_back(point_7);
    }
};

class Triangle : public Polygon {
public:
    Point o_center() {
      std::vector<Point> d = getVertices();
      Line line_1(d[0], d[1]);
      Line line_2(d[1], d[2]);
      Point point_1 = mid(d[0], d[1]);
      Line normal_1 = normal(line_1, point_1);
      Point point_2 = mid(d[1], d[2]);
      Line normal_2 = normal(line_2, point_2);
      Point center = crossing(normal_1, normal_2);
      return center;
    }

    Circle circumscribedCircle() {
      std::vector<Point> d = getVertices();
      Point center = o_center();
      double a = dist(d[0], d[1]);
      double b = dist(d[1], d[2]);
      double c = dist(d[2], d[0]);
      double p = (a + b + c) / 2;
      double r = (a * b * c) / (4 * sqrt(p * (p - a) * (p - b) * (p - c)));
      Circle ans(center, r);
      return ans;
    }

    Point incenter() {
      std::vector<Point> d = getVertices();
      Line line_1(d[0], d[1]);
      Line line_2(d[1], d[2]);
      Line line_3(d[2], d[0]);
      Line bisectrix_1 = bisectrix(line_1, line_2);
      Line bisextrix_2 = bisectrix(line_2, line_3);
      Point center = crossing(bisectrix_1, bisextrix_2);
      return center;
    }


    Circle inscribedCircle() {
      std::vector<Point> d = getVertices();
      Point center = incenter();
      double a = dist(d[0], d[1]);
      double b = dist(d[1], d[2]);
      double c = dist(d[2], d[0]);
      double p = (a + b + c) / 2;
      double r = sqrt((p - a) * (p - b) * (p - c) / p);
      Circle ans(center, r);
      return ans;
    }

    Point centroid() {
      std::vector<Point> d = getVertices();
      Point mid_1 = mid(d[0], d[1]);
      Point mid_2 = mid(d[1], d[2]);
      Line line_1(mid_1, d[2]);
      Line line_2(mid_2, d[0]);
      return crossing(line_1, line_2);
    }

    Point orthocenter() {
      std::vector<Point> d = getVertices();
      Line line_1(d[0], d[1]);
      Line line_2(d[1], d[2]);
      Line normal_1 = normal(line_1, d[2]);
      Line normal_2 = normal(line_2, d[0]);
      return crossing(normal_1, normal_2);
    }

    Line EulerLine() {
      Line ans(orthocenter(), incenter());
      return ans;
    }

    Circle ninePointsCircle() {
      Point point_1 = orthocenter();
      Point point_2 = o_center();
      double r = circumscribedCircle().radius();
      return Circle(mid(point_1, point_2), r / 2);
    }

    Triangle(const Point &point_1, const Point &point_2, const Point &point_3) {
      vertices_.push_back(point_1);
      vertices_.push_back(point_2);
      vertices_.push_back(point_3);
    }
};

