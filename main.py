import math as m


class Point_GeographicalCoordinates:
    def __init__(self, latitudine, longitudine):
        """
        :param latitudine: latitude of a given point
        :param longitudine: longitude of a given point
            --> A point is defined by 2 values (coordinates) : latitude & longitude;
        """

        self.latitudine = latitudine
        self.longitudine = longitudine


class SphericalCoordinate:

    def __init__(self, gc1: tuple, gc2: tuple):
        """
        :param gc1: is a tuple that contains the first Point (coordinates: latitude + longitude)
        :param gc2: is a tuple that contains the second Point (coordinates: latitude + longitude)
        """
        self.gc1 = gc1
        self.gc2 = gc2

    @property  # @property ---> transform a method into an argument
    def azimuth(self):
        """
        :return: A (in degrees)
            A = the azimtuh between 2 points
                * A = atan2(sin(ec_1), cos(ec_2)) --- atan2 is much more precisely than atan
                * A is initial in radians
        """

        try:
            l = self.gc2[1] - self.gc1[1]
            ec_1 = m.cos(m.radians(self.gc2[0])) * m.sin(m.radians(l))
            ec_2 = m.sin(m.radians(self.gc2[0])) * m.cos(m.radians(self.gc1[0])) - \
                   m.cos(m.radians(self.gc2[0])) * m.sin(m.radians(self.gc1[0])) * m.cos(m.radians(l))
            A = m.degrees(m.atan2(ec_1, ec_2))

            return A
        except:
            return "ERROR: incorrect input for azimuth "

    @property
    def distance(self):
        """
        :return:   the distance (in meters) between the two locations
        """

        # a = 6378388
        #
        # c = 0.999
        #
        # phi_0 = 45.9
        # lambda_0 = 25.39213
        #
        #
        # e_2 = 0.0066943800
        #
        #
        # X_point1 = (a * m.cos(m.radians(self.gc1[0])) * m.cos(m.radians(self.gc1[1]))) / (m.sqrt(1 - e_2 * (m.sin(m.radians(self.gc1[0])))**2))
        # Y_point1 = (a * m.cos(m.radians(self.gc1[0])) * m.sin(m.radians(self.gc1[1]))) / (m.sqrt(1 - e_2 * (m.sin(m.radians(self.gc1[0])))**2))
        #
        # X_point2 = (a * m.cos(m.radians(self.gc2[0])) * m.cos(m.radians(self.gc2[1]))) / (m.sqrt(1 - e_2 * (m.sin(m.radians(self.gc2[0])))**2))
        # Y_point2 = (a * m.cos(m.radians(self.gc2[0])) * m.sin(m.radians(self.gc2[1]))) / (m.sqrt(1 - e_2 * (m.sin(m.radians(self.gc2[0])))**2))
        #
        # return ((m.sqrt((X_point2-X_point1)**2 + (Y_point2-Y_point1)**2)))
        #
        # l1 = self.gc1[1] - lambda_0
        # l2 = self.gc2[1] - lambda_0
        #
        # x1 = (m.sin(m.radians(self.gc1[0])) * m.cos(m.radians(phi_0)) -
        #       m.cos(m.radians(self.gc1[0])) * m.sin(m.radians(phi_0)) * m.cos(m.radians(l1))) \
        #      / (1+m.sin(m.radians(self.gc1[0])) * m.sin(m.radians(phi_0)) +
        #         m.cos(m.radians(self.gc1[0])) * m.cos(m.radians(phi_0)) * m.cos(m.radians(l1)))
        #
        # y1 = (m.cos(m.radians(self.gc1[0]))*m.sin(m.radians(l1))) / (1+m.sin(m.radians(self.gc1[0])) * m.sin(m.radians(phi_0)) +
        #         m.cos(m.radians(self.gc1[0])) * m.cos(m.radians(phi_0)) * m.cos(m.radians(l1)))
        #
        #
        # x2 =  (m.sin(m.radians(self.gc2[0])) * m.cos(m.radians(phi_0)) -
        #       m.cos(m.radians(self.gc2[0])) * m.sin(m.radians(phi_0)) * m.cos(m.radians(l2))) \
        #      / (1+m.sin(m.radians(self.gc2[0])) * m.sin(m.radians(phi_0)) +
        #         m.cos(m.radians(self.gc2[0])) * m.cos(m.radians(phi_0)) * m.cos(m.radians(l2)))
        #
        # y2 = (m.cos(m.radians(self.gc2[0]))*m.sin(m.radians(l2))) / (1+m.sin(m.radians(self.gc2[0])) * m.sin(m.radians(phi_0)) +
        #         m.cos(m.radians(self.gc2[0])) * m.cos(m.radians(phi_0)) * m.cos(m.radians(l2)))
        #
        # X2 = x2*c + 500000
        # X1 = x1*c + 500000
        #
        # Y1 = y1*c + 500000
        # Y2 = y2*c + 500000
        #
        # dx = X2-X1
        # dy = Y2-Y1
        # #
        # #
        # return (m.sqrt(4*(a**2)*(dx**2 + dy**2)))

        try:
            R = 6378137

            dlon = m.radians(self.gc2[1]) - m.radians(self.gc1[1])
            dlat = m.radians(self.gc2[0]) - m.radians(self.gc1[0])

            a = m.sin(dlat / 2) ** 2 + m.cos(self.gc1[0]) * m.cos(self.gc2[0]) * m.sin(dlon / 2) ** 2
            c = 2 * m.atan(m.sqrt(a) / m.sqrt(1 - a))

            return R * c
        except:
            return "ERROR: incorrect input for distance"

    def endPosition(self, gc: tuple, azimuth, arc_distance):
        """
        :param gc:  Geographical Coordinates for the start point
        :param azimuth: azimuth value in radians
        :param arc_distance: distnace between start point and end point

        :return: the location on a large arc of a circle with the given:
                            initial location, azimuth and distance of the arc (in radians)
        """
        try:
            R = 6378137
            A = m.radians(azimuth)
            d = arc_distance  # (1000)

            a = 6378137
            b = 6356752.3142

            Z = 0.2 * m.atan(m.sqrt(a - b) / m.sqrt(a + b))

            end_latitude = m.asin(m.cos(A) * m.tan(Z) * m.cos(gc[0]) * m.cos(gc[1]))
            end_latitude = round(end_latitude, 4)

            end_longitude = m.asin((m.sin(Z) * m.sin(A)) / (m.cos(0)))
            end_longitude = round(end_longitude, 3)

            return (end_latitude, end_longitude)
        except:
            return "ERROR: incorrect input for 'endPosision' method "

    def interpolate(self, amount):
        """
        :param amount:  interpolation factor

                --> The interpolation factor represents a value in the range [0,1]. If it has the value 0 then
                       the starting position is returned and if it has the value 1 the end position is returned.

        :return: -	Returns an interpolated value between the two locations (coordinates) along a large arc.
        """

        try:
            delta_latitude = self.gc2[0] - self.gc1[0]
            delta_longitude = self.gc2[1] - self.gc1[1]

            return (amount * delta_latitude, amount * delta_longitude)
        except:
            return "ERROR: incorrect input for 'interpolate' method"


# Here i created the poins :
A = Point_GeographicalCoordinates(0.0000, 0.0000)
B = Point_GeographicalCoordinates(0.0000, 45.0000)

# Here i created the Spherical Coordinates for points:
SphericalCoordinate_of_points = SphericalCoordinate((A.latitudine, A.longitudine), (B.latitudine, B.longitudine))

"""  Results:  """
print('\n')
print("* Exercitiul 1: \n    Output : ", SphericalCoordinate_of_points.azimuth)

print("* Exercitiul 2: \n    Output : ", SphericalCoordinate_of_points.distance)

print("* Exercitiul 3: \n    Output : ",
      SphericalCoordinate_of_points.endPosition((0.000, 0.000), SphericalCoordinate_of_points.azimuth, 1000) )

print("* Exercitiul 4: \n    Output : ", SphericalCoordinate_of_points.interpolate(0.5) )